! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module L2GPData                 ! Creation, manipulation and I/O for L2GP Data
!=============================================================================
  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use DUMP_0, only: DUMP
  use Hdf, only: DFACC_READ, DFACC_CREATE, DFACC_RDWR, &
    & DFNT_CHAR8, DFNT_FLOAT32, DFNT_INT32, DFNT_FLOAT64
  use HDFEOS, only: SWATTACH, SWDETACH, SWINQDIMS
  use Intrinsic ! "units" type literals, beginning with L_
  use MLSCommon, only: I4, R4, R8, DEFAULTUNDEFINEDVALUE
  use MLSFiles, only: FILENOTFOUND, &
    & HDFVERSION_4, HDFVERSION_5, WILDCARDHDFVERSION, WRONGHDFVERSION, &
    & MLS_HDF_VERSION, MLS_INQSWATH, MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF, &
    & MLS_EXISTS
  use MLSHDFEOS, only: mls_swattach, mls_swdetach
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning, MLSMSG_Debug
  use MLSStrings, only: Capitalize,  ints2Strings, lowercase, &
    & strings2Ints
  use MLSStringLists, only: ExtractSubString, &
    & GetStringHashElement, GetStringElement, GetUniqueList, &
    & list2array, NumStringElements, &
    & StringElementNum
  use OUTPUT_M, only: BLANKS, OUTPUT
  use PCFHdr, only: GA_VALUE_LENGTH, GlobalAttributes_T, GlobalAttributes, &
    & he5_readglobalattr, he5_writeglobalattr
  use STRING_TABLE, only: DISPLAY_STRING

  implicit none

  private
  public :: L2GPData_T
  public :: L2GPNameLen
  public :: AddL2GPToDatabase, AppendL2GPData, cpL2GPData, &
    & DestroyL2GPContents,  DestroyL2GPDatabase, diff, Dump, &
    & ExpandL2GPDataInPlace, &
    & ReadL2GPData, SetupNewL2GPRecord, WriteL2GPData

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       & "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character(len=*), parameter, private :: ModuleName = &
       & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  interface DIFF
    module procedure DiffL2GPData
    module procedure DiffL2GPFiles
  end interface

  interface DIFFSTATS
    module procedure DiffStatsInt
    !module procedure DiffStatsSingle
    !module procedure DiffStatsDouble
  end interface

  interface DUMP
    module procedure DUMP_L2GP
    module procedure DUMP_L2GP_DataBase
    module procedure DumpL2GP_attributes_hdf5
  end interface

  interface ReadL2GPData
    module procedure ReadL2GPData_fileID
    module procedure ReadL2GPData_fileName
  end interface

  interface AppendL2GPData
    module procedure AppendL2GPData_fileID
    module procedure AppendL2GPData_fileName
  end interface

  interface cpL2GPData
    module procedure cpL2GPData_fileID
    module procedure cpL2GPData_fileName
  end interface

  ! This module defines datatypes and gives basic routines for storing and
  ! manipulating L2GP data.
  ! It is prepared to handle io for both file versions: hdfeos2 and hdfeos5

!     c o n t e n t s
!     - - - - - - - -

! L2GPData_T              The l2gp data type; holds all data for one swath

! AddL2GPToDatabase       Adds an l2gp data type to a database
! DestroyL2GPContents     Deallocates all the arrays allocated for an L2GP
! DestroyL2GPDatabase     Destroys an L2GP database
! Dump                    Reveals info about an L2GP or a database of L2GP
!                            to any desired level of detail
! ExpandL2GPDataInPlace   Adds more profiles to an existing L2GP
! AppendL2GPData          Appends L2GP onto end of existing swath file
! ReadL2GPData            Reads L2GP from existing swath file
! SetupNewL2GPRecord      Allocates arrays for a new L2GP
! WriteL2GPData           Writes an L2GP into a swath file

  ! First some local parameters
  ! Assume L2GP files w/o explicit hdfVersion field are this
  ! 4 corresponds to hdf4, 5 to hdf5 in L2GP, L2AUX, etc. 
  integer, parameter :: L2GPDEFAULT_HDFVERSION = HDFVERSION_5

  ! r4 corresponds to sing. prec. :: same as stored in files
  integer, public, parameter :: rgp = r4

  ! How long may the list of swath names grow (~80 x max num. of swaths/file)
  integer, public, parameter :: MAXNUMSWATHPERFILE = 250
  integer, public, parameter :: MAXSWATHNAMESBUFSIZE = 80*MAXNUMSWATHPERFILE

  ! TRUE means we can avoid using unlimited dimension and its time penalty
  logical, public            :: AVOIDUNLIMITEDDIMS = .true.

  integer, parameter :: CHARATTRLEN = 255   ! was GA_VALUE_LENGTH
  real, parameter    :: UNDEFINED_VALUE = DEFAULTUNDEFINEDVALUE !-999.99 ! Same as %template%badvalue
  integer, parameter :: L2GPNameLen = 80
  integer, parameter :: NumGeolocFields = 10
  integer, parameter :: MAXNLEVELS = 1000

   ! The following are the current data fields
   character (len=*), parameter :: DATA_FIELD1 = 'L2gpValue'
   character (len=*), parameter :: DATA_FIELD2 = 'L2gpPrecision'
   character (len=*), parameter :: DATA_FIELD3 = 'Status'
   character (len=*), parameter :: DATA_FIELD4 = 'Quality'
   ! This is the above except for Status which has proved troublesome
   character (len=*), parameter :: DATA_FIELDS = &
     & 'L2gpValue,L2gpPrecision,Quality'

   ! The following are the current geolocation fields
   character (len=*), parameter :: GEO_FIELD1 = 'Latitude'
   character (len=*), parameter :: GEO_FIELD2 = 'Longitude'
   character (len=*), parameter :: GEO_FIELD3 = 'Time'
   character (len=*), parameter :: GEO_FIELD4 = 'LocalSolarTime'
   character (len=*), parameter :: GEO_FIELD5 = 'SolarZenithAngle'
   character (len=*), parameter :: GEO_FIELD6 = 'LineOfSightAngle'
   character (len=*), parameter :: GEO_FIELD7 = 'OrbitGeodeticAngle'
   character (len=*), parameter :: GEO_FIELD8 = 'ChunkNumber'
   character (len=*), parameter :: GEO_FIELD9 = 'Pressure'
   character (len=*), parameter :: GEO_FIELD10= 'Frequency'
   ! These are the above except for 8,9,10 which have special attributes
   character (len=*), parameter :: GEO_FIELDS = &
     & 'Latitude,Longitude,LocalSolarTime,SolarZenithAngle,LineOfSightAngle' // &
     & ',OrbitGeodeticAngle'

   ! The following are the dimension names according to the mls spec
   character (len=*), parameter :: DIM_NAME1 = 'nTimes'
   character (len=*), parameter :: DIM_NAME2 = 'nLevels'
   character (len=*), parameter :: DIM_NAME3 = 'nFreqs'
   ! An alternate name for DIM_NAME3 used by HIRDLS
   character (len=*), parameter :: HRD_DIM_NAME3 = 'nChans'
   character (len=*), parameter :: DIM_NAME12 = 'nLevels,nTimes'
   character (len=*), parameter :: DIM_NAME123 = 'nFreqs,nLevels,nTimes'
   ! These are for the new max_dimlist parameter added to  SWdefgfld.
   ! this one is for non-extendible dimensions
   character (len=*), parameter :: MAX_DIML = ' '
   character (len=*), parameter :: UNLIM = 'Unlim'
   ! This is for cases where the time dimension is extendible
   character (len=*), parameter :: MAX_DIML1 = UNLIM
   character (len=*), parameter :: MAX_DIML12 = 'nLevels,Unlim'
   character (len=*), parameter :: MAX_DIML123 = 'nFreqs,nLevels,Unlim'

!   INTEGER,PARAMETER::CHUNKFREQS=13,CHUNKLEVELS=17,CHUNKTIMES=9,CHUNK4=1
   
   character (len=*), parameter :: DEFAULTMAXDIM = UNLIM
   integer, parameter :: DEFAULT_CHUNKRANK = 1
   integer, dimension(7), parameter :: DEFAULT_CHUNKDIMS = (/ 1,1,1,1,1,1,1 /)
   integer, parameter :: HDFE_AUTOMERGE = 1     ! MERGE FIELDS WITH SHARE DIM
   integer, parameter :: HDFE_NOMERGE = 0       ! don't merge
  ! The following, if true, says to encode strings as ints
  ! before swapi write; also decode ints to strings after read
  ! otherwise, try swapi read, write directly with strings
  ! logical, parameter :: USEINTS4STRINGS = .false.  
  
  ! So far, the nameIndex component of the main data type is never set
  logical, parameter :: NAMEINDEXEVERSET = .false.  
  
  ! Do you want to write file and swath attributes when you append values
  logical, parameter :: APPENDSWRITEATTRIBUTES = .true.  

  ! Do you want to pre-fill arrays with MissingValue by default?
  logical, parameter :: ALWAYSFILLWITHMISSINGVALUE = .true.  

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

     character (LEN=L2GPNameLen) :: name ! Typically the swath name.
     integer :: nameIndex       ! Used by the parser to keep track of the data
     integer :: QUANTITYTYPE = 0   ! E.g., l_temperature; (Where is this used?)

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

     ! We always write this data field
     ! However, because some l2gp files were created incorrectly,
     ! It takes a special forcing option to read it
     ! character (len=1), pointer, dimension(:) :: status=>NULL()
     ! Now we've changed our minds: status will be a 4-byte integer
     integer(i4), pointer, dimension(:) :: status=>NULL()
     !                (status is a reserved word in F90)
     real (rgp), pointer, dimension(:) :: quality=>NULL()
     ! Both the above dimensioned (nTimes)

     ! The dimensions for the quantity (if, e.g., coming from l2cf)
     ! character(len=CHARATTRLEN)        :: DIM_Names = '' ! ','-separated
     ! character(len=CHARATTRLEN)        :: DIM_Units = '' ! ','-separated
     ! character(len=CHARATTRLEN)        :: VALUE_Units = '' 
     ! These are the fill/missing values for all arrays except status
     real (rgp)                        :: MissingValue = UNDEFINED_VALUE
     integer(i4)                       :: MissingStatus = 513 ! 512 + 1
    ! Vertical coordinate
    character(len=8) :: verticalCoordinate ! E.g. 'Pressure', or 'Theta'
    ! integer :: verticalCoordinate ! The vertical coordinate used.  These
                                  ! are l_lits of the type t_VGridCoord
                                  ! defined in Init_Tables_Module.
  end type L2GPData_T

  ! Print debugging stuff?
  logical, parameter :: DEEBUG = .false.  
  logical, parameter ::SWATHLEVELMISSINGVALUE = .false. ! Make it swath attr?
  logical, parameter ::READINGSTATUSBYDEFAULT = .true.  ! Change if bombing

contains ! =====     Public Procedures     =============================

  !------------------------------------------  SetupNewL2GPRecord  -----
  subroutine SetupNewL2GPRecord ( l2gp, nFreqs, nLevels, nTimes, nTimesTotal, &
    & FillIn)

    ! This routine sets up the arrays for an l2gp datatype.

    ! Dummy arguments
    type (L2GPData_T), intent(inout)  :: l2gp
    integer, intent(in), optional :: nFreqs            ! Dimensions
    integer, intent(in), optional :: nLevels           ! Dimensions
    integer, intent(in), optional :: nTimes            ! Dimensions
    integer, intent(in), optional :: nTimesTotal       ! Dimensions
    logical, intent(in), optional :: FillIn    ! Fill with MissingValue

    ! Local variables
    integer :: useNFreqs, useNLevels, useNTimes, useNTimesTotal
    logical :: myFillIn

    if ( present(nFreqs) ) then
       useNFreqs=nFreqs
    else
       useNFreqs=0
    end if

    if ( present(nLevels) ) then
       useNLevels=nLevels
    else
       useNLevels=0
    end if

    if ( present(nTimes) ) then
       useNTimes=nTimes
    else
       useNTimes=0              ! Default to empty l2gp
    end if

    if ( present(nTimesTotal) ) then
       useNTimesTotal=nTimesTotal
    else
       useNTimesTotal=0              ! Default to empty l2gp
    end if

    if ( present(FillIn) ) then
      myFillIn = FillIn
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
    if ( .not. myFillIn ) return
    l2gp%pressures = l2gp%MissingValue
    l2gp%frequency = l2gp%MissingValue
    l2gp%latitude = l2gp%MissingValue
    l2gp%longitude = l2gp%MissingValue
    l2gp%solarTime = l2gp%MissingValue
    l2gp%losAngle = l2gp%MissingValue
    l2gp%geodAngle = l2gp%MissingValue
    l2gp%time = l2gp%MissingValue
    l2gp%chunkNumber = l2gp%MissingValue
    l2gp%l2gpValue = l2gp%MissingValue
    l2gp%l2gpPrecision = l2gp%MissingValue
    l2gp%status = l2gp%MissingStatus ! l2gp%MissingValue
    ! l2gp%status = ' '
    l2gp%quality = l2gp%MissingValue

  end subroutine SetupNewL2GPRecord

  !-----------------------------------------  DestroyL2GPContents  -----
  subroutine DestroyL2GPContents ( L2GP )

    ! This routine deallocates all the arrays allocated above.

    ! Dummy arguments
    type (L2GPData_T), intent(inout) :: L2GP

    ! Executable code

    call deallocate_test ( l2gp%pressures,         "l2gp%pressures",         ModuleName )
    call deallocate_test ( l2gp%latitude,          "l2gp%latitude",          ModuleName )
    call deallocate_test ( l2gp%longitude,         "l2gp%longitude",         ModuleName )
    call deallocate_test ( l2gp%solarTime,         "l2gp%solarTime",         ModuleName )
    call deallocate_test ( l2gp%solarZenith,       "l2gp%solarZenith",       ModuleName )
    call deallocate_test ( l2gp%losAngle,          "l2gp%losAngle",          ModuleName )
    call deallocate_test ( l2gp%geodAngle,         "l2gp%geodAngle",         ModuleName )
    call deallocate_test ( l2gp%chunkNumber,       "l2gp%chunkNumber",       ModuleName )
    call deallocate_test ( l2gp%time,              "l2gp%time",              ModuleName )
    call deallocate_test ( l2gp%frequency,         "l2gp%frequency",         ModuleName )
    call deallocate_test ( l2gp%l2gpValue,         "l2gp%l2gpValue",         ModuleName )
    call deallocate_test ( l2gp%l2gpPrecision,     "l2gp%l2gpPrecision",     ModuleName )
    call deallocate_test ( l2gp%status,            "l2gp%status",            ModuleName )
    call deallocate_test ( l2gp%quality,           "l2gp%quality",           ModuleName )
    l2gp%nTimes = 0
    l2gp%nTimesTotal = 0
    l2gp%nLevels = 0
    l2gp%nFreqs = 0

  end subroutine DestroyL2GPContents

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
      & l2gp%l2gpPrecision, l2gp%status, l2gp%quality )
    
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
    
    ! Deallocate the old arrays
    call DestroyL2GPContents(templ2gp)

  end subroutine ExpandL2GPDataInPlace

  !-------------------------------------------  AddL2GPToDatabase  -----
  integer function AddL2GPToDatabase( DATABASE, ITEM )

    ! This function adds an l2gp data type to a database of said types,
    ! creating a new database if it doesn't exist.  The result value is
    ! the size -- where L2gp is put.

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

    ! Dummy argument
    type (l2GPData_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer :: l2gpIndex, status

    if ( associated(database) ) then
       do l2gpIndex = 1, SIZE(database)
          call DestroyL2GPContents ( database(l2gpIndex) )
       end do
       deallocate ( database, stat=status )
       if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & MLSMSG_deallocate // "database" )
    end if
  end subroutine DestroyL2GPDatabase

  ! ---------------------- ReadL2GPData_fileID  -----------------------------

  subroutine ReadL2GPData_fileID(L2FileHandle, swathname, l2gp, numProfs, &
       firstProf, lastProf, hdfVersion, HMOT, ReadStatus)
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
    logical, optional, intent(in) :: ReadStatus

    ! Local
    integer :: myhdfVersion
    character :: my_hmot
    
    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L2GPDEFAULT_HDFVERSION
    endif
    my_hmot = 'M'
    if ( present(hmot)) my_hmot = Capitalize(hmot)
    ! Check for valid hmot
    select case (my_hmot)
    case ('H')
    case ('M')
    case ('O')
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Unable to Read OMI L2GPData" )
    case ('T')
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Unable to Read TES L2GPData" )
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Unrecognized instrument key passed to ReadL2GPData: "// my_hmot )
    end select
    !print*,"In readl2gpdata: first/last prof=",firstProf, lastProf
    if (myhdfVersion == HDFVERSION_4) then
      call ReadL2GPData_hdf(L2FileHandle, swathname, l2gp, my_hmot, myhdfVersion,&
        & numProfs, firstProf, lastProf, ReadStatus)
    elseif (myhdfVersion /= HDFVERSION_5) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Unrecognized hdfVersion passed to ReadL2GPData" )
    else
      call ReadL2GPData_hdf(L2FileHandle, swathname, l2gp, my_hmot, myhdfVersion,&
        & numProfs, firstProf, lastProf, ReadStatus)
    endif
    !print*,"In readl2gpdata: first/last/read prof=",firstProf,&
    !  lastProf,numProfs
  end subroutine ReadL2GPData_fileID

  ! ---------------------- ReadL2GPData_fileName  -----------------------------

  subroutine ReadL2GPData_fileName(fileName, swathname, l2gp, numProfs, &
       firstProf, lastProf, hdfVersion, HMOT, ReadStatus)
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
    logical, optional, intent(in) :: ReadStatus

    ! Local
    integer :: L2FileHandle
    integer :: record_length
    integer :: status
    integer :: the_hdfVersion
    
    ! Executable code
    the_hdfVersion = mls_hdf_version(FileName, hdfVersion)
    if ( the_hdfVersion == FILENOTFOUND ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'File not found; make sure the name and path are correct' &
        & // trim(fileName) )

    L2FileHandle = mls_io_gen_openF('swopen', .TRUE., status, &
       & record_length, DFACC_READ, FileName=FileName, &
       & hdfVersion=hdfVersion, debugOption=.false. )
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to open L2gp file: " // trim(FileName) // ' for reading')
    call ReadL2GPData_fileID(L2FileHandle, swathname, l2gp, numProfs=numProfs, &
       & firstProf=firstProf, lastProf=lastProf, hdfVersion=the_hdfVersion, &
       & hmot=hmot)
    status = mls_io_gen_closeF('swclose', L2FileHandle, FileName=FileName, &
      & hdfVersion=hdfVersion)
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2gp file: " // trim(FileName) // ' after reading')
  end subroutine ReadL2GPData_fileName

  ! ------------------- ReadL2GPData_hdf ----------------

  subroutine ReadL2GPData_hdf(L2FileHandle, swathname, l2gp, HMOT, hdfVersion,&
    & numProfs, firstProf, lastProf, ReadStatus)
  use HDFEOS5, only: HE5_swattach, HE5_swdetach, HE5_SWINQDIMS, HE5_swfldinfo
  use MLSHDFEOS, only: mls_swdiminfo, mls_swrdfld
    !------------------------------------------------------------------------

    ! This routine reads an L2GP file, returning a filled data structure and the !
    ! number of profiles read.

    ! All the ReadStatus harrumphing is because
    ! (1) An earlier version always core dumped
    ! (2) Even this version dumps core while reading l2gp files which stored
    !     the 'Status' field as 32-bit floats instead of chars
    ! Therefore, unless you supply the optional arg ReadStatus=.true.,
    ! it will skip reading the troublesome datafield
    ! Naturally, you would do that only for correctly-formatted l2gp files
    ! Arguments

    character (LEN=*), intent(IN) :: swathname ! Name of swath
    integer, intent(IN) :: L2FileHandle ! Returned by HE5_swopen
    integer, intent(IN) :: hdfVersion ! Returned by HE5_swopen
    integer, intent(IN), optional :: firstProf, lastProf ! Defaults to first and last
    type( L2GPData_T ), intent(OUT) :: l2gp ! Result
    character, intent(in) :: HMOT   ! 'H' 'M'(def) 'O' 'T'
    integer, intent(OUT),optional :: numProfs ! Number actually read
    logical, optional, intent(in) :: ReadStatus

    ! Local Parameters
    character (LEN=*), parameter :: SZ_ERR = 'Failed to get size of &
         &dimension '
    character (LEN=*), parameter :: MLSMSG_INPUT = 'Error in input argument '
    character (LEN=*), parameter :: MLSMSG_L2GPRead = 'Unable to read L2GP &
                                                     &field:'
    ! Local Variables
    character (len=80) :: DF_Name
    character (len=80) :: DF_Precision
    character (LEN=80) :: list
    character (LEN=80) :: dimlist
    character (LEN=80) :: maxdimlist
    character (LEN=8)  :: maxdimName
    integer :: rank
    integer, dimension(7) :: numberType
    integer, dimension(7) :: flddims
    character (LEN=480) :: msr

    integer :: alloc_err, first, freq, lev, nDims, size, swid, status
    integer :: start(3), stride(3), edge(3), dims(3)
    integer :: nFreqs, nLevels, nTimes, nFreqsOr1, nLevelsOr1, myNumProfs
    logical :: firstCheck, lastCheck, timeIsUnlim

    real(r4), pointer, dimension(:) :: REALFREQ
    real(r4), pointer, dimension(:) :: REALSURF
    real(r4), pointer, dimension(:) :: REALPROF
    real(r4), pointer, dimension(:,:,:) :: REAL3
    logical :: dontfail
    logical :: ReadingStatus
    logical :: deeBugHere
!  How to deal with status and columnTypes? swrfld fails
!  with char data on Linux with HDF4. With HDF5 we may or may not need to
!  Have recourse to ints2Strings and strings2Ints if USEINTS4STRINGS
!    character (LEN=8), allocatable :: the_status_buffer(:)
!    character (LEN=L2GPNameLen), allocatable :: the_status_buffer(:)
!    integer, allocatable, dimension(:,:) :: string_buffer

    deeBugHere = DEEBUG     ! .or. .true.
    nullify ( realFreq, realSurf, realProf, real3 )
    
    ! Don't fail when trying to read an mls-specific field 
    ! if the file is from another Aura instrument
    dontfail = (HMOT /= 'M')
    ReadingStatus = READINGSTATUSBYDEFAULT ! was .false.
    if ( present(ReadStatus) ) ReadingStatus = ReadStatus
    ! Attach to the swath for reading
    ! print*," in ReadL2GPData_hdf: first/last=",firstprof,lastprof
    l2gp%Name = swathname
    
    ! print *, 'Trying to read he5_swattach to read'
    select case (HMOT)
    case ('H')
      swid = mls_SWattach(L2FileHandle, 'HIRDLS', hdfVersion=hdfVersion)
      DF_Name = TRIM(l2gp%Name)
      DF_Precision = TRIM(l2gp%Name) // 'Precision'
      l2gp%MissingValue = -999.  ! This is a HIRDLS-specific setting
    case ('M')
      swid = mls_SWattach(L2FileHandle, l2gp%Name, hdfVersion=hdfVersion)
      DF_Name = DATA_FIELD1
      DF_Precision = DATA_FIELD2
    case default
    end select
    if (swid == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &attach to hdfeos2/5 swath interface for reading' // trim(swathname))

    ! Get dimension information

    lev = 0
    freq = 0

    if( hdfVersion == HDFVERSION_4 ) then
      nDims = swinqdims(swid, list, dims)
    else
      nDims = HE5_SWinqdims(swid, list, dims)
    endif
    if ( deeBugHere ) print *, 'HMOT: ', HMOT
    if ( deeBugHere ) print *, 'swathName: ', l2gp%name
    if ( deeBugHere ) print *, 'dimlist: ', trim(list)
    if ( deeBugHere ) print *, 'ndims: ', ndims
    if ( deeBugHere ) print *, 'dims: ', dims
    if (nDims == -1) call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'Failed to get dimension information on hdfeos5 swath ' // trim(swathname))
    if ( index(list,'nLevels') /= 0 ) lev = 1
    if ( index(list,'Freq') /= 0 ) freq = 1
!    if ( index(list,'Unlim') /= 0 ) then 
!      timeIsUnlim = .TRUE.
!    else
!      timeIsUnlim = .FALSE.
!    endif

    size = mls_swdiminfo(swid, 'nTimes', hdfVersion=hdfVersion)
    if ( hdfVersion == HDFVERSION_5 ) then
      ! This will be wrong if timeIsUnlim .eq. .TRUE. . 
      ! HE5_SWdiminfo returns 1 instead of the right answer.
      ! print *,"just called mls_swdiminfo(nTimes) : ", size
      status = HE5_swfldinfo(swid, trim(DF_Name), rank, flddims, &
        & numberType, dimlist, maxdimlist)
      ! print *, 'rank: ', rank
      ! print *, 'flddims: ', flddims(1:rank)
      ! print *, 'numberType: ', numberType(1:rank)
      ! print *, 'dimlist: ', dimlist
      ! print *, 'maxdimlist: ', maxdimlist
      call GetStringHashElement (dimlist, &
       & maxdimlist, 'nTimes', &
       & maxDimName, .false.)
   !   print *, 'maxdimName: ', maxdimName
      if ( maxDimName == 'Unlim' ) then
        status = StringElementNum(dimlist, 'nTimes', .false.)
       ! print *, 'elem num: ', status
        if ( status > 0 .and. status <= rank ) size = max(size, flddims(status))
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
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
       else
          first = firstProf
       endif

    else

       first = 0

    endif

    if (lastCheck) then

       if (lastProf < first) then
          msr = MLSMSG_INPUT // 'lastProf'
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
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
    ! l2gp%status = ' ' ! So it has a value.
    if ( ReadingStatus) &
      & status = mls_swrdfld( swid, 'Status',start(3:3),stride(3:3),edge(3:3),&
      & l2gp%status, hdfVersion=hdfVersion, dontfail=.true. )

    if ( HMOT == 'M' ) then
      status = mls_SWrdfld(swid, 'Quality', start(3:3), stride(3:3),&
        edge(3:3),realProf, hdfVersion=hdfVersion, dontfail=dontfail)
      l2gp%quality = realProf
    endif

    ! Deallocate local variables
    call Deallocate_test ( realProf, 'realProf', ModuleName )
    call Deallocate_test ( realSurf, 'realSurf', ModuleName )
    call Deallocate_test ( realFreq, 'realFreq', ModuleName )
    call Deallocate_test ( real3, 'real3', ModuleName )

    !  After reading, detach from HE5_SWath interface
    status = mls_SWdetach(swid, hdfVersion=hdfVersion)
    if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &detach from swath interface after reading.')
    !print*," leaving ReadL2GPData_hdf: first/last/read=",&
    !  firstprof,lastprof,myNumProfs
    ! Set numProfs if wanted
    if (present(numProfs)) numProfs=myNumProfs


    !-----------------------------
  end subroutine ReadL2GPData_hdf
  !-----------------------------

  ! --------------------------------------  OutputL2GP_createFile_hdf  -----
  subroutine OutputL2GP_createFile_hdf (l2gp, L2FileHandle, hdfVersion, &
    & swathName, fileName, nLevels, notUnlimited, compressTimes)

  use HDFEOS5, only: HE5_SWdetach, &
    & HE5S_UNLIMITED_F, &
    & HE5T_NATIVE_CHAR, HE5T_NATIVE_DOUBLE, HE5T_NATIVE_INT, HE5T_NATIVE_FLOAT
  use MLSHDFEOS, ONLY : mls_swcreate, mls_dfldsetup, mls_gfldsetup, &
    & mls_swdefdim
    ! Brief description of subroutine
    ! This subroutine sets up the structural definitions in an empty L2GP file.

    ! Arguments

    integer, intent(in) :: L2FileHandle ! From swopen
    type( L2GPData_T ), intent(inout) :: l2gp
    integer, intent(in) :: hdfVersion
    character (LEN=*), optional, intent(IN) :: swathName ! Defaults to l2gp%swathName
    character (LEN=*), optional, intent(IN) :: fileName
    integer, optional, intent(in) :: nLevels
    logical, optional, intent(in) :: notUnlimited   !               as nTimes
    logical, optional, intent(in) :: compressTimes  ! don't store nTimesTotal
    ! Parameters

    character (len=*), parameter :: DIM_ERR = 'Failed to define dimension '
    character (len=*), parameter :: GEO_ERR = &
         & 'Failed to define geolocation field '
    character (len=*), parameter :: DAT_ERR = 'Failed to define data field '

    ! Variables

    character (len=480) :: MSR
    character (len=132) :: NAME   ! From l2gp%name
    character (len=32) :: MYDIM1, MYDIM12, MYDIM123

    ! THESE ARE HDF5 CHUNKS, _NOT_ MLS ALONG-TRACK PROCESSING CHUNKS 
    integer,dimension(7)::CHUNK_DIMS
    integer::CHUNK_RANK
    integer::CHUNKTIMES,CHUNKFREQS,CHUNKLEVELS

    integer :: SWID, STATUS
    logical :: myNotUnlimited
    logical :: mycompressTimes
    integer, external :: he5_SWgetfill
    real :: fillValue
    ! integer, external :: he5_swsetfill

    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    myNotUnlimited = .false.
    if ( present ( notUnlimited ) ) myNotUnlimited = notUnlimited
    mycompressTimes = .false.
    if ( present (compressTimes ) ) mycompressTimes = compressTimes
    
    ! Work out the chunking
    if ( myNotUnlimited ) then
      chunkTimes = max ( min ( 120, l2gp%nTimes ), 1 )
    else
      chunktimes = 120      ! was 1
    endif
    chunkfreqs = max ( l2gp%nFreqs, 1)
    if(present(nLevels))then
       chunklevels = nLevels
    else
      chunklevels = min(l2gp%nLevels, 500)     ! was .., 5)
      chunklevels = max(chunklevels, 1)
    endif
    
    ! Create the swath within the file
    ! print*,"Creating swath called ",name

    if ( DEEBUG )  print *, 'About to sw_create ', TRIM(name)
    if ( present(filename) .and. DEEBUG ) print *, 'file name ', TRIM(filename)
    swid = mls_SWcreate(L2FileHandle, trim(name), &
      & filename=filename, hdfVersion=hdfVersion)
    !print*,"Swath ",name,"has SW id :",swid
    if ( swid == -1 ) then
       msr = 'Failed to create swath ' // TRIM(name) &
        & // ' (maybe has the same name as another swath in this file?)'
    end if

    ! Define dimensions

    if ( hdfVersion == HDFVERSION_5 .and. .not. myNotUnlimited ) then
      ! Defining special "unlimited dimension called UNLIM
      ! print*,"Defined Unlim with size", HE5S_UNLIMITED_f
      ! status = HE5_SWdefdim(swid, UNLIM, HE5S_UNLIMITED_F)
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
    if ( DEEBUG ) then
      print *, 'myDim1 ', myDim1
      print *, 'myDim12 ', myDim12
      print *, 'myDim123 ', myDim123
      print *, 'nTimes ', l2gp%nTimes
      print *, 'nTimesTotal ', l2gp%nTimesTotal
      print *, 'nLevels ', l2gp%nLevels
      print *, 'nFreqs ', l2gp%nFreqs
    endif
    if ( mycompressTimes ) then
      status = mls_swdefdim(swid, 'nTimes', max(l2gp%nTimes,1), &
        & hdfVersion=hdfVersion)
    else
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
      & hdfVersion=hdfVersion, iFill=int(l2gp%MissingValue))

    if ( l2gp%nLevels > 0 ) then

      status = mls_gfldsetup(swid, 'Pressure', 'nLevels', MAX_DIML, &
      & DFNT_FLOAT32, HDFE_NOMERGE, 0, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)
    end if

    if ( l2gp%nFreqs > 0 ) then

      status = mls_gfldsetup(swid, 'Frequency', 'nFreqs', MAX_DIML, &
      & DFNT_FLOAT32, HDFE_NOMERGE, 0, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)
    end if

    ! Define data fields using above dimensions

    if ( (l2gp%nFreqs > 0) .and. (l2gp%nLevels > 0) ) then
       chunk_rank=3
       chunk_dims(1:3)=(/ CHUNKFREQS,CHUNKLEVELS,CHUNKTIMES /)

      status = mls_dfldsetup(swid, 'L2gpValue', 'nFreqs,nLevels,nTimes', &
      & MYDIM123, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

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
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

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
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

      status = mls_dfldsetup(swid, 'L2gpPrecision', 'nTimes', &
      & MYDIM1, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

    end if

    chunk_rank=1
    chunk_dims(1)=CHUNKTIMES
    status = mls_dfldsetup(swid, 'Status', 'nTimes', &
    & MYDIM1, &
    & DFNT_INT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
    & hdfVersion=hdfVersion, iFill=l2gp%MissingStatus)
    ! & hdfVersion=hdfVersion, iFill=int(l2gp%MissingValue))
    ! & DFNT_CHAR8, HDFE_NOMERGE, chunk_rank, chunk_dims, &

    chunk_rank=1
    chunk_dims(1)=CHUNKTIMES
    status = mls_dfldsetup(swid, 'Quality', 'nTimes', &
    & MYDIM1, &
    & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
    & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

    ! Detach from the HE5_SWath interface.This stores the swath info within the
    ! file and must be done before writing or reading data to or from the
    ! swath. (May be un-necessary for HDF5 -- test program works OK without.)
    ! 
    status = mls_SWdetach(swid, hdfVersion=hdfVersion)
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Failed to detach from swath interface after definition.' )
    end if

    !--------------------------------------
  end subroutine OutputL2GP_createFile_hdf
  !--------------------------------------

  !-----------------------------------------  OutputL2GP_writeGeo_hdf  -----
  subroutine OutputL2GP_writeGeo_hdf (l2gp, l2FileHandle, hdfVersion, &
    & swathName,offset)

  use HDFEOS5, only: HE5_swattach, HE5_swdetach
  use MLSHDFEOS, only: mls_swwrfld
    ! Brief description of subroutine
    ! This subroutine writes the geolocation fields to an L2GP output file.

    ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp
    integer, intent(in) :: l2FileHandle ! From swopen
    integer, intent(in) :: hdfVersion
    character (len=*), intent(IN), optional :: swathName ! Defaults->l2gp%name
    integer,intent(IN),optional::offset
    ! Parameters

    character (len=*), parameter :: WR_ERR = &
         & 'Failed to write geolocation field '
    
    ! Variables

    character (len=132) :: name ! Either swathName or l2gp%name
    
    integer :: status, swid,myOffset
    integer :: start(2), stride(2), edge(2)

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

    ! print *, 'Trying to he5_swattach to write geo'
    swid = mls_SWattach (l2FileHandle, name, hdfVersion=hdfVersion)

    ! Write data to the fields

    stride = 1
    start = myOffset ! Please do not set to zero
    edge(1) = l2gp%nTimes
    ! print *, 'Writing geolocation fields'
    ! print *, 'start', start
    ! print *, 'stride', stride
    ! print *, 'edge', edge
    ! print *, 'shape(Latitude)', shape(l2gp%latitude)

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

    status = mls_SWwrfld(swid, 'ChunkNumber', start, stride, edge, &
         l2gp%chunkNumber, hdfVersion=hdfVersion)

    if ( l2gp%nLevels > 0 ) then
       edge(1) = l2gp%nLevels
       start(1)=0 ! needed because offset may have made this /=0
       status = mls_SWwrfld(swid, 'Pressure', start, stride, edge, &
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
            & 'Failed to detach from swath interface' )
    end if

    !------------------------------------
  end subroutine OutputL2GP_writeGeo_hdf
  !------------------------------------

  !----------------------------------------  OutputL2GP_writeData_hdf  -----
  subroutine OutputL2GP_writeData_hdf(l2gp, l2FileHandle, hdfVersion, &
    & swathName,offset)

  use HDFEOS5, only: HE5_swattach, HE5_swdetach
  use MLSHDFEOS, only: mls_swwrfld
    ! Brief description of subroutine
    ! This subroutine writes the data fields to an L2GP output file.
    ! For now, you have to write all of l2gp, but you can choose to write
    ! it at some offset into the file
    ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp
    integer, intent(in) :: l2FileHandle ! From swopen
    integer, intent(in) :: hdfVersion
    character (len=*), intent(IN), optional :: swathName ! Defaults->l2gp%name
    integer,intent(IN),optional::offset
    ! Parameters

    character (len=*), parameter :: WR_ERR = 'Failed to write data field '

    ! Variables

    character (len=480) :: msr
    character (len=132) :: name     ! Either swathName or l2gp%name

    integer :: status,myOffset
    integer :: start(3), stride(3), edge(3)
    integer :: swid
    
!  How to deal with status and columnTypes? swrfld fails
!  with char data on Linux
!  Have recourse to ints2Strings and strings2Ints if USEINTS4STRINGS
!    character (LEN=8), allocatable :: the_status_buffer(:)
!    character (LEN=L2GPNameLen), allocatable :: the_status_buffer(:)
!    integer, allocatable, dimension(:,:) :: string_buffer

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

    start = 0
    stride = 1
    start(3)= myOffset ! Please do not set to zero
    edge(1) = l2gp%nFreqs
    edge(2) = l2gp%nLevels
    edge(3) = l2gp%nTimes
    ! print *, 'Trying to he5_swattach to write data'
    swid = mls_SWattach (l2FileHandle, name, hdfVersion=hdfVersion)
    if ( l2gp%nFreqs > 0 ) then
       ! Value and Precision are 3-D fields
       if (DEEBUG) print *, 'start, stride, edge ', start, stride, edge
       status = mls_SWwrfld(swid, 'L2gpValue', start, stride, edge, &
            & reshape(real(l2gp%l2gpValue), (/size(l2gp%l2gpValue)/)), &
            & hdfVersion=hdfVersion )
       status = mls_SWwrfld(swid, 'L2gpPrecision', start, stride, edge, &
            & reshape(real(l2gp%l2gpPrecision), (/size(l2gp%l2gpPrecision)/)), &
            & hdfVersion=hdfVersion )

    else if ( l2gp%nLevels > 0 ) then
       ! Value and Precision are 2-D fields
      
       if (DEEBUG) print *, 'start, stride, edge: ', start(2:3), stride(2:3), edge(2:3)
       status = mls_SWwrfld( swid, 'L2gpValue', start(2:3), stride(2:3), &
            edge(2:3), real(l2gp%l2gpValue(1,:,:) ), hdfVersion=hdfVersion)
       status = mls_SWwrfld( swid, 'L2gpPrecision', start(2:3), stride(2:3), &
            edge(2:3), real(l2gp%l2gpPrecision(1,:,:) ), hdfVersion=hdfVersion)
    else

       ! Value and Precision are 1-D fields
       if (DEEBUG) print *, 'start, stride, edge: ', start(3), stride(3), edge(3)
       status = mls_SWwrfld( swid, 'L2gpValue', start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpValue(1,1,:) ), hdfVersion=hdfVersion)
       status = mls_SWwrfld( swid, 'L2gpPrecision', start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpPrecision(1,1,:) ), hdfVersion=hdfVersion)
    end if

    ! 1-D status & quality fields

    status = mls_swwrfld(swid, 'Status', start(3:3), stride(3:3), edge(3:3), &
       &   l2gp%status, hdfVersion=hdfVersion, dontfail=.true.)

    status = mls_SWwrfld(swid, 'Quality', start(3:3), stride(3:3), edge(3:3), &
         real(l2gp%quality), hdfVersion=hdfVersion)

    !     Detach from the swath interface.

    status = mls_SWdetach(swid, hdfVersion=hdfVersion)
    if(DEEBUG) print *, 'Detached from swid ', swid
    if(DEEBUG) print *, 'file handle ', l2FileHandle
    if(DEEBUG) print *, 'status ', status
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach from swath interface' )
    end if


  !-------------------------------------
  end subroutine OutputL2GP_writeData_hdf
  !-------------------------------------

  !----------------------------------------  OutputL2GP_attributes_hdf5  -----
  subroutine OutputL2GP_attributes_hdf5(l2gp, l2FileHandle, swathName)

  use HDFEOS5, only: HE5T_NATIVE_INT, HE5T_NATIVE_REAL, HE5T_NATIVE_DOUBLE, &
    & HE5_SWattach, HE5_SWdetach, MLS_charType
  use he5_swapi, only: he5_swwrattr, he5_swwrlattr
  use MLSHDFEOS, only: mls_swwrattr, mls_swwrlattr
  use PCFHdr, only:  he5_writeglobalattr
    ! Brief description of subroutine
    ! This subroutine writes the attributes for an l2gp
    ! These include
    ! Global attributes which are file level
    ! Swath level
    ! Geolocation field attributes
    ! Data field attributes
    ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp
    integer, intent(in) :: l2FileHandle ! From swopen
    character (len=*), intent(IN), optional :: swathName ! Defaults->l2gp%name
    ! Parameters

    character (len=*), parameter :: WR_ERR = 'Failed to write attribute field '

    ! Variables

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
    integer :: rgp_type
    integer :: status
    integer :: swid
    character(len=CHARATTRLEN), dimension(NumGeolocFields) :: theTitles
    character(len=CHARATTRLEN), dimension(NumGeolocFields) :: theUnits
    character(len=CHARATTRLEN) :: field_name
    character(len=CHARATTRLEN) :: species_name
    character(len=CHARATTRLEN) :: units_name
    character(len=CHARATTRLEN) :: abbr_uniq_fdef
    character(len=CHARATTRLEN) :: expnd_uniq_fdef
    
    ! Begin
    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    if ( rgp == r4 ) then
      rgp_type = HE5T_NATIVE_REAL
    elseif ( rgp == r8 ) then
      rgp_type = HE5T_NATIVE_DOUBLE
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & 'Attributes have unrecognized numeric data type; should be r4 or r8')
    endif
    call List2Array(GeolocationTitles, theTitles, .true.)
    call List2Array(GeolocationUnits, theUnits, .true.)

    ! - -   G l o b a l   A t t r i b u t e s   - -
    call he5_writeglobalattr(l2FileHandle)

    swid = mls_SWattach (l2FileHandle, name, hdfVersion=HDFVERSION_5)
    
    !   - -   S w a t h   A t t r i b u t e s   - -
    status = he5_swwrattr(swid, trim(l2gp%verticalCoordinate), &
      & rgp_type, size(l2gp%pressures), &
      & l2gp%pressures)
    field_name = l2gp%verticalCoordinate ! 'Pressure'
    status = mls_swwrattr(swid, 'VerticalCoordinate', MLS_CHARTYPE, 1, &
      & field_name)
    if ( SWATHLEVELMISSINGVALUE ) &
      & status = he5_swwrattr(swid, 'MissingValue', rgp_type, 1, &
      & (/ real(l2gp%MissingValue, rgp) /) )
    
    !   - -   G e o l o c a t i o n   A t t r i b u t e s   - -
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
        call GetStringHashElement (GeolocationTitles, &
          & GeoUniqueFieldDefinition, trim(theTitles(field)), &
          & abbr_uniq_fdef, .false.)
        call GetStringHashElement (UniqueFieldDefKeys, &
          & UniqueFieldDefValues, trim(abbr_uniq_fdef), &
          & expnd_uniq_fdef, .false.)
        status = mls_swwrlattr(swid, trim(theTitles(field)), 'Title', &
          & MLS_CHARTYPE, 1, theTitles(field))
        status = mls_swwrlattr(swid, trim(theTitles(field)), 'Units', &
          & MLS_CHARTYPE, 1, theUnits(field))

        if ( trim(theTitles(field)) == 'Time' ) then
          status = he5_swwrlattr(swid, trim(theTitles(field)), &
            & 'MissingValue', &
            & HE5T_NATIVE_DOUBLE, 1, (/ real(l2gp%MissingValue, r8) /) )
        elseif ( trim(theTitles(field)) == 'ChunkNumber' ) then
          status = he5_swwrlattr(swid, trim(theTitles(field)), &
            & 'MissingValue', &
            & HE5T_NATIVE_INT, 1, (/ int(l2gp%MissingValue) /) )
        else
          status = he5_swwrlattr(swid, trim(theTitles(field)), &
            & 'MissingValue', &
            & rgp_type, 1, (/ real(l2gp%MissingValue, rgp) /) )
        endif
        status = mls_swwrlattr(swid, trim(theTitles(field)), &
          & 'UniqueFieldDefinition', &
          & MLS_CHARTYPE, 1, trim(expnd_uniq_fdef))
      endif
    enddo
    !   - -   D a t a   A t t r i b u t e s   - -
    ! call GetQuantityAttributes ( l2gp%quantityType, &
    !  & units_name, expnd_uniq_fdef)
    field_name = Name
    species_name = name
    isColumnAmt = ( index(species_name, 'Column') > 0 )
    if ( isColumnAmt ) then
      call ExtractSubString(Name, species_name, 'Column', 'wmo')
    endif
    call GetStringHashElement (lowercase(Species), &
      & SpUniqueFieldDefinition, trim(lowercase(species_name)), &
      & abbr_uniq_fdef, .false.)
    call GetStringHashElement (UniqueFieldDefKeys, &
      & UniqueFieldDefValues, trim(abbr_uniq_fdef), &
      & expnd_uniq_fdef, .false.)
    if ( expnd_uniq_fdef == '' .or. expnd_uniq_fdef == ',' ) &
      & expnd_uniq_fdef = 'MLS-Specific'
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
    status = mls_swwrlattr(swid, 'L2gpValue', 'Title', &
      & MLS_CHARTYPE, 1, field_name)
    status = mls_swwrlattr(swid, 'L2gpValue', 'Units', &
      & MLS_CHARTYPE, 1, units_name)
    status = he5_swwrlattr(swid, 'L2gpValue', 'MissingValue', &
      & rgp_type, 1, (/ real(l2gp%MissingValue, rgp) /) )
    status = mls_swwrlattr(swid, 'L2gpValue', &
      & 'UniqueFieldDefinition', &
      & MLS_CHARTYPE, 1, trim(expnd_uniq_fdef))
    status = mls_swwrlattr(swid, 'L2gpPrecision', 'Title', &
      & MLS_CHARTYPE, 1, trim(field_name)//'Precision')
    status = mls_swwrlattr(swid, 'L2gpPrecision', 'Units', &
      & MLS_CHARTYPE, 1, units_name)
    status = he5_swwrlattr(swid, 'L2gpPrecision', 'MissingValue', &
      & rgp_type, 1, (/ real(l2gp%MissingValue, rgp) /) )
    status = mls_swwrlattr(swid, 'L2gpPrecision', &
      & 'UniqueFieldDefinition', &
      & MLS_CHARTYPE, 1, trim(expnd_uniq_fdef))

    ! ('Status' data field newly written)
    status = mls_swwrlattr(swid, 'Status', 'Title', &
      & MLS_CHARTYPE, 1, trim(field_name)//'Status')
    status = mls_swwrlattr(swid, 'Status', 'Units', &
      & MLS_CHARTYPE, 1, 'NoUnits')
    status = he5_swwrlattr(swid, 'Status', 'MissingValue', &
      & HE5T_NATIVE_INT, 1, (/ l2gp%MissingStatus /) )
      ! & HE5T_NATIVE_INT, 1, (/ int(l2gp%MissingValue) /) )

    ! status = mls_swwrlattr(swid, 'Status', 'MissingValue', &
    !   & MLS_CHARTYPE, 1, ' ' )
    status = mls_swwrlattr(swid, 'Status', &
      & 'UniqueFieldDefinition', &
      & MLS_CHARTYPE, 1, 'MLS-Specific')
    
    status = mls_swwrlattr(swid, 'Quality', 'Title', &
      & MLS_CHARTYPE, 1, trim(field_name)//'Quality')
    status = mls_swwrlattr(swid, 'Quality', 'Units', &
      & MLS_CHARTYPE, 1, units_name)
    status = he5_swwrlattr(swid, 'Quality', 'MissingValue', &
      & rgp_type, 1, (/ real(l2gp%MissingValue, rgp) /) )
    status = mls_swwrlattr(swid, 'Quality', &
      & 'UniqueFieldDefinition', &
      & MLS_CHARTYPE, 1, 'MLS-Specific')
    
    status = mls_SWdetach(swid, hdfVersion=HDFVERSION_5)
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach  from swath interface' )
    end if


  !-------------------------------------
  end subroutine OutputL2GP_attributes_hdf5
  !-------------------------------------

  !----------------------------------------  SetL2GP_aliases  -----
  subroutine SetL2GP_aliases(l2gp, l2FileHandle, swathName)

  use HDFEOS5, only: HE5_SWATTACH, HE5_SWSETALIAS, HE5_SWDETACH
  use SDPToolkit, only: PGS_S_SUCCESS
    ! Arguments
    integer, intent(IN) :: l2FileHandle ! From swopen
    type (L2GPData_T), intent(INOUT) :: l2gp
    character (LEN=*), optional, intent(IN) ::swathName!default->l2gp%swathName
    
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
    ! print *, 'Trying to he5_swattach to set alias'
    sw_id = mls_swattach(l2FileHandle, trim(name), hdfVersion=HDFVERSION_5)
    if ( sw_id < 1 ) then 
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & "Error in attaching swath for setting alias." )
    end if
    returnStatus = he5_SWsetalias(sw_id, TYPE2FIELDNAME, trim(name))
    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & "Error in setting alias from " // TYPE2FIELDNAME // &
          & ' to ' // trim(name) )
    end if
    returnStatus = he5_SWsetalias(sw_id, TYPE2PRECISIONNAME, &
     & trim(name) // 'Precision')
    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & "Error in setting alias from " // TYPE2PRECISIONNAME // &
          & ' to ' // trim(name) // 'Precision' )
    end if
    returnStatus = mls_SWdetach(sw_id, hdfVersion=HDFVERSION_5)
    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & "Error in detaching swath for setting alias." )
    end if
  !-------------------------------------
  end subroutine SetL2GP_aliases
  !-------------------------------------


  ! This subroutine is an amalgamation of the last three
  ! Should be renamed CreateAndWriteL2GPData
  subroutine WriteL2GPData(l2gp, l2FileHandle, swathName, filename, hdfVersion, &
    & notUnlimited)

    ! Arguments

    integer, intent(IN) :: l2FileHandle ! From swopen
    type (L2GPData_T), intent(INOUT) :: l2gp
    character (LEN=*), optional, intent(IN) ::swathName!default->l2gp%swathName
    character (LEN=*), optional, intent(IN) ::fileName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: notUnlimited
    ! Exectuable code

    ! Local
    integer :: myhdfVersion

    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L2GPDEFAULT_HDFVERSION
    endif

    call OutputL2GP_createFile_hdf (l2gp, l2FileHandle, myhdfVersion, &
      & swathName, filename, notUnlimited=notUnlimited)
    call OutputL2GP_writeGeo_hdf (l2gp, l2FileHandle, myhdfVersion, &
      & swathName)
    call OutputL2GP_writeData_hdf (l2gp, l2FileHandle, myhdfVersion, &
      & swathName)
    if (myhdfVersion == HDFVERSION_5) then
      if ( DEEBUG ) print *, 'Outputting attributes'
      call OutputL2GP_attributes_hdf5 (l2gp, l2FileHandle, swathName)
      if ( DEEBUG ) print *, 'Setting aliases'
      call SetL2GP_aliases (l2gp, l2FileHandle, swathName)
    endif

  end subroutine WriteL2GPData
  !-------------------------------------------------------------

  subroutine AppendL2GPData_fileID(l2gp, l2FileHandle, &
    & swathName, filename, offset, lastProfile, TotNumProfs, hdfVersion, &
    & createSwath)
    ! sticks l2gp into the swath swathName in the file pointed at by
    ! l2FileHandle,starting at the profile number "offset" (First profile
    ! in the file has offset==0). If this runs off the end of the swath, 
    ! it is lengthened automagically. 
    ! This call has been altered recently, so that it can be used to create
    ! a swath as well as adding to one. 

    ! Arguments

    integer, intent(IN) :: l2FileHandle ! From swopen

    ! This is a L2GPData_T structure containing all the data to be written
    type (L2GPData_T), intent(INOUT) :: l2gp
    ! This is the name the swath is given in the file. By default it is
    ! the name contained in l2gp
    character (LEN=*), optional, intent(IN) ::swathName!default->l2gp%swathName
    character (LEN=*), optional, intent(IN) ::fileName
    ! This (offset) is the point in the swath at which the data is written. 
    ! First profile in the file has offset==0. If the swath in the file is 
    ! shorter than offset + ( num of profiles in l2gp) then it grows by magic
    integer, intent(IN), optional::offset
    ! TotNumProfs is a new argument. It seems only to be used if we are 
    ! creating a swath, rather than adding to one. In that case I guess
    ! it is the total number of profiles in the swath created. I also 
    ! guess that this is done so that we can avoid growing and re-growing 
    ! the swath.
    integer, intent(IN), optional::TotNumProfs
    integer, intent(IN), optional::lastProfile
    integer, optional, intent(in) :: hdfVersion ! better be 4 or 5!
    logical, intent(in), optional :: createSwath
    ! Local
    integer :: actual_ntimes
    integer :: myhdfVersion
    integer :: status
    logical :: swath_exists
    integer :: swathid
    integer :: myLastProfile
    character (len=L2GPNameLen) :: myswathName

    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L2GPDEFAULT_HDFVERSION
    endif
    ! Optional args should _ONLY_ appear like this, inside an IF 
    ! block that checks if they are present. HCP replaced one occurance
    ! of TotNumProfs with myLastProfile as they are equal if TotNumProfs
    ! is present and cause a crash if it is not. Another occurrance seemed 
    ! wrong anyway, so I commented it out.
    if (present(lastProfile)) then
      myLastProfile = lastProfile 
    elseif (present(TotNumProfs)) then
      myLastProfile = TotNumProfs 
    else
      myLastProfile = L2GP%nTimesTotal
    endif
    myswathName = l2gp%name
    if ( present(swathName) ) myswathName = swathName
    
    if ( present(createSwath) ) then
      swath_exists = .not. createSwath
      ! print *, 'createSwath: ', createSwath
    else
      ! print *, 'Uh-oh, calling mls_swattach'
      swathid = mls_swattach(L2FileHandle, trim(myswathName), &
        & hdfVersion=myhdfVersion, DONTFAIL=.true.)
      swath_exists = ( swathid > 0 )
      if ( swath_exists ) then
        status = mls_swdetach(swathid, hdfVersion=myhdfVersion)
        if ( status /= 0 ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, & 
          & 'Failed to detach from swath in AppendL2GPData_fileID')
      endif
    endif

    if ( swath_exists ) then
      if(DEEBUG) print *, 'OK, swath already exists'
    else
      ! Must create swath in file w/o disturbing other swaths
      if(DEEBUG) print *, 'Must create swath'
      if(DEEBUG) print *, 'Will have ', myLastProfile, ' profiles'
      if(DEEBUG) print *, 'instead of ', l2gp%nTimes, ' profiles'
      actual_ntimes = l2gp%nTimes
      ! if ( present(TotNumProfs) ) l2gp%nTimes = TotNumProfs
      select case (myhdfVersion)
      case (HDFVERSION_4)
        ! Currently force unlimited, remove the .false. .and. to allow limited
        if ( present(TotNumProfs) ) l2gp%nTimes = TotNumProfs
        call OutputL2GP_createFile_hdf (l2gp, L2FileHandle, myhdfVersion, &
          & myswathName, filename, notUnlimited=.false. .and. present(totNumProfs))
        l2gp%nTimes = actual_ntimes
      case (HDFVERSION_5)
        ! By default allow limited; 
        ! may force unlimited by setting avoidUnlimitedDims to FALSE
        ! if ( present(TotNumProfs) ) l2gp%nTimes = TotNumProfs
        call OutputL2GP_createFile_hdf (l2gp, L2FileHandle, myhdfVersion, &
          & myswathName, filename,&
          & notUnlimited=(avoidUnlimitedDims .and. present(totNumProfs)) )
          ! & myswathName, filename, notUnlimited=present(totNumProfs))
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName, &
         & 'Illegal hdf version in AppendL2GPData_fileName')
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
    ! This line caused an error because TotNumProfs is optional and 
    ! is not here checked for present-ness. 
    ! if ( offset == TotNumProfs .or. size(l2gp%l2gpValue,3) == 0 ) then
    if ( size(l2gp%l2gpValue,3) == 0 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & "No profiles in this chunk" )

    else
      ! actual_ntimes = l2gp%nTimes
      ! l2gp%nTimes = max(myLastProfile - offset + 1, 1)
      call OutputL2GP_writeGeo_hdf (l2gp, l2FileHandle, myHDFVersion, &
        & myswathName, offset)
      call OutputL2GP_writeData_hdf (l2gp, l2FileHandle, myHDFVersion, &
        & myswathName, offset)
      select case ( myHDFVersion )
      case ( HDFVERSION_4 )
      case ( HDFVERSION_5 )
        if ( .not. swath_exists .and. APPENDSWRITEATTRIBUTES) then
          call OutputL2GP_attributes_hdf5 (l2gp, l2FileHandle, swathName)
          call SetL2GP_aliases (l2gp, l2FileHandle, swathName)
        end if
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Unrecognized hdfVersion passed to AppendL2GPData" )
      end select
      ! l2gp%nTimes = actual_ntimes
    end if

  end subroutine AppendL2GPData_fileID

  ! ---------------------- AppendL2GPData_fileName  ---------------------------

  subroutine AppendL2GPData_fileName(l2gp, fileName, &
    & swathname, offset, lastProfile, TotNumProfs, hdfVersion)
    !------------------------------------------------------------------------

    ! Given a file name,
    ! This routine does an append operation (see AppendL2GPData_fileID)
    ! If the file doesn't exist yet, hopefully, it'll create it

    ! Arguments

    character (len=*), intent(in) :: fileName ! Name of swath
    type( l2GPData_T ), intent(inout) :: l2gp ! Result
    character (LEN=*), optional, intent(IN) ::swathName!default->l2gp%swathName
    integer,intent(IN),optional :: offset
    integer, intent(IN), optional::lastProfile
    integer,intent(IN),optional :: TotNumProfs
    integer, optional, intent(in) :: hdfVersion
    ! logical, optional, intent(in) :: clean   ! Not implemented yet

    ! Local
    integer :: L2FileHandle
    integer :: record_length
    integer :: status
    integer :: the_hdfVersion
    logical :: myClean
    logical :: file_exists
    integer :: file_access
    integer :: actual_ntimes
    
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
    L2FileHandle = mls_io_gen_openF('swopen', .TRUE., status, &
       & record_length, file_access, FileName=FileName, &
       & hdfVersion=the_hdfVersion, debugOption=.false. )
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to open L2gp file: " // trim(FileName) // ' for appending')
    call AppendL2GPData_fileID(l2gp, L2FileHandle, swathname, &
      & FileName, offset, lastProfile=lastProfile, totNumProfs=totNumProfs, &
      & hdfVersion=the_hdfVersion)
    status = mls_io_gen_closeF('swclose', L2FileHandle, FileName=FileName, &
      & hdfVersion=the_hdfVersion)
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2gp file: " // trim(FileName) // ' after appending')
  end subroutine AppendL2GPData_fileName

  ! ---------------------- cpL2GPData_fileID  ---------------------------

  subroutine cpL2GPData_fileID(file1, file2, swathList, &
    & hdfVersion1, hdfVersion2, notUnlimited, swathList2, ReadStatus)
    !------------------------------------------------------------------------

    ! Given file names file1 and file2,
    ! This routine copies swathList from 1 to 2

    ! Arguments

    integer, intent(in)           :: file1 ! handle of file 1
    integer, intent(in)           :: file2 ! handle of file 2
    character (len=*), intent(in) :: swathList
    integer, intent(in)           :: hdfVersion1
    integer, intent(in)           :: hdfVersion2
    logical, optional, intent(in) :: notUnlimited
    logical, optional, intent(in) :: ReadStatus
    character (len=*), optional, intent(in) :: swathList2

    ! Local variables
    logical, parameter            :: countEmpty = .true.
    type (L2GPData_T) :: l2gp
    integer :: i
    integer :: noSwaths
    character (len=L2GPNameLen) :: swath
    character (len=L2GPNameLen) :: swath2
    
    ! Executable code
    noSwaths = NumStringElements(trim(swathList), countEmpty)
    if ( noSwaths < 1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No swaths cp to file--unable to count swaths in ' // trim(swathList) )
    endif
    if ( present(swathList2) ) &
      & noSwaths = min(noSwaths, NumStringElements(trim(swathList2), countEmpty))
    if ( noSwaths < 1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No swaths cp to file--unable to count swaths in ' // trim(swathList2) )
    endif
    ! Loop over swaths in file 1
    do i = 1, noSwaths
      call GetStringElement (trim(swathList), swath, i, countEmpty )
      if ( present(swathList2) ) then
        call GetStringElement (trim(swathList2), swath2, i, countEmpty )
      else
        swath2 = swath
      endif
      ! Allocate and fill l2gp
      if ( DEEBUG ) print *, 'Reading swath from file: ', trim(swath)
      call ReadL2GPData ( file1, trim(swath), l2gp, &
           & hdfVersion=hdfVersion1, ReadStatus=ReadStatus )
      if ( DEEBUG ) then
        print *, 'Writing swath to file: ', trim(swath)
        print *, 'l2gp%nFreqs:  ', l2gp%nFreqs
        print *, 'l2gp%nLevels: ', l2gp%nLevels
        print *, 'l2gp%nTimes:  ', l2gp%nTimes
        print *, 'shape(l2gp%l2gpvalue):  ', shape(l2gp%l2gpvalue)
      endif
      ! Write the filled l2gp to file2
      call WriteL2GPData(l2gp, file2, trim(swath2), hdfVersion=hdfVersion2, &
        & notUnlimited=notUnlimited)
      call DestroyL2GPContents ( l2gp )
    enddo
       
  end subroutine cpL2GPData_fileID

  ! ---------------------- cpL2GPData_fileName  ---------------------------

  subroutine cpL2GPData_fileName(file1, file2, &
    & create2, hdfVersion1,  hdfVersion2, swathList, swathList2, &
    & notUnlimited, andGlAttributes, ReadStatus)
    !------------------------------------------------------------------------

    ! Given file names file1 and file2,
    ! This routine copies all the l2gpdata from 1 to 2
    ! (see cpL2GPData_fileID)
    ! If file2 doesn't exist yet, or if create2 is TRUE, it'll create it

    ! Arguments

    character (len=*), intent(in) :: file1 ! Name of file 1
    character (len=*), intent(in) :: file2 ! Name of file 2
    logical, optional, intent(in) :: create2
    logical, optional, intent(in) :: andGlAttributes
    integer, optional, intent(in) :: hdfVersion1
    integer, optional, intent(in) :: hdfVersion2
    character (len=*), optional, intent(in) :: swathList
    character (len=*), optional, intent(in) :: swathList2
    logical, optional, intent(in) :: notUnlimited
    logical, optional, intent(in) :: ReadStatus

    ! Local
    integer :: File1Handle
    integer :: File2Handle
    integer :: record_length
    integer :: status
    integer :: the_hdfVersion1
    integer :: the_hdfVersion2
    logical :: file_exists
    integer :: file_access
    integer :: listsize
    integer :: noSwaths
    character (len=MAXSWATHNAMESBUFSIZE) :: mySwathList
    type(GlobalAttributes_T) :: gAttributes
    type(GlobalAttributes_T) :: gAttributesOriginal
    character(len=40)        :: ProcessLevel
    integer                  :: DayofYear
    double precision         :: TAI93At0zOfGranule
    logical                  :: myandGlAttributes
    
    ! Executable code
    the_hdfVersion1 = L2GPDEFAULT_HDFVERSION
    if ( present(hdfVersion1) ) the_hdfVersion1 = hdfVersion1
    file_exists = ( mls_exists(trim(File1)) == 0 )
    if ( .not. file_exists ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
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
    if ( present(swathList) ) then
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
    File1Handle = mls_io_gen_openF('swopen', .TRUE., status, &
       & record_length, DFACC_READ, FileName=File1, &
       & hdfVersion=the_hdfVersion1, debugOption=.false. )
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to open L2gp file: " // trim(File1) // ' for cp-ing')

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
    File2Handle = mls_io_gen_openF('swopen', .TRUE., status, &
       & record_length, file_access, FileName=File2, &
       & hdfVersion=the_hdfVersion2, debugOption=.false. )
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to open L2gp file: " // trim(File2) // ' for cping')
    if ( DEEBUG ) then
      print *, 'About to cp from file1 to file2: ', File1Handle, File2Handle
      print *, trim(mySwathList)
    endif
    ! Maybe copy global attributes, too
    if ( myandGlAttributes ) then
      call he5_readglobalattr (File1Handle, gAttributes, &
        & ProcessLevel, DayofYear, TAI93At0zOfGranule, status)
      if ( status == 0 ) then
        ! Unfortunately, he5_writeglobalattr writes class-level data
        ! so must store original version so can copy new data into it 
        gAttributesOriginal = GlobalAttributes
        GlobalAttributes = gAttributes
        call he5_writeglobalattr (File2Handle)
        ! Before leaving must copy original data back
        GlobalAttributes = gAttributesOriginal
      endif
    endif
    call cpL2GPData_fileID(File1Handle, File2Handle, &
      & mySwathList, the_hdfVersion1, the_hdfVersion2, &
      & notUnlimited=notUnlimited, swathList2=swathList2, &
      & ReadStatus=ReadStatus )
    if ( DEEBUG ) print *, 'About to close File1Handle: ', File1Handle
    status = mls_io_gen_closeF('swclose', File1Handle, FileName=File1, &
      & hdfVersion=the_hdfVersion1, debugOption=.false.)
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2gp file: " // trim(File1) // ' after cping')
    if ( DEEBUG ) print *, 'About to close File2Handle: ', File2Handle
    ! status = mls_io_gen_closeF('swclose', File2Handle, FileName=File2, &
    !  & hdfVersion=the_hdfVersion, debugOption=.true.)
    status = mls_io_gen_closeF('swclose', File2Handle, &
      & hdfVersion=the_hdfVersion2, debugOption=.false.)
    if ( status /= 0 ) then
      print *, 'status returned from mls_io_gen_closeF: ', status
      print *, 'WRONGHDFVERSION: ', WRONGHDFVERSION
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2gp file: " // trim(File2) // ' after cping')
    endif
  end subroutine cpL2GPData_fileName

  ! ------------------------------------------ DiffL2GPData ------------
  subroutine DiffL2GPData ( L2gp1, L2gp2, &
    & Details, wholeArray, stats, rms, ignoreBadChunks )
    ! Show diff between l2gp1 and l2gp2 down to level of Details
    ! Assumes fields of each already allocated
    ! (If not, then why are you trying to show differences?)
    ! (Couldn't you at least print warning and return if not?)
    
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
    logical, intent(in), optional :: STATS   ! if TRUE, just print stats
    logical, intent(in), optional :: RMS     ! if TRUE, just print mean, rms
    logical, intent(in), optional :: WHOLEARRAY   ! if TRUE, print anyway
    logical, intent(in), optional :: IGNOREBADCHUNKS   ! if TRUE, ignore
                                                     ! instances where geod bad
    ! Local variables
    integer :: ierr
    integer :: MYDETAILS
    real(r8) :: FillValue
    integer :: ChunkFillValue
    logical :: ShapesDontMatch
    logical :: badChunks
    integer :: instance
    integer :: badInstances
    ! Executable code
    myDetails = 1
    if ( present(details) ) myDetails = details
    ShapesDontMatch = .false.
    badChunks = .false.
    badInstances = 0
      !FillValue = real(l2gp1%MissingValue, r8)
      !ChunkFillValue = int(l2gp1%MissingValue)

      if ( trim(l2gp1%name) /= trim(l2gp2%name) ) then
        call output('(1) name: ' // trim(l2gp1%name), advance='yes')
        call output('(2) name: ' // trim(l2gp2%name), advance='yes')
      endif
      if ( myDetails < -1 ) return
      if ( L2gp1%MissingValue /= L2gp2%MissingValue ) then
        call output(' (1) MissingValue = ', advance='no')
        call output(L2gp1%nTimes, advance='yes')
        call output(' (2) MissingValue = ', advance='no')
        call output(L2gp2%nTimes, advance='yes')
      endif
      if ( L2gp1%nTimes /= L2gp2%nTimes ) then
        call output(' (1) nTimes = ', advance='no')
        call output(L2gp1%nTimes, advance='yes')
        call output(' (2) nTimes = ', advance='no')
        call output(L2gp2%nTimes, advance='yes')
        ShapesDontMatch = .true.
      endif
      if ( L2gp1%nLevels /= L2gp2%nLevels ) then
        call output(' (1) nLevels = ', advance='no')
        call output(L2gp1%nLevels, advance='yes')
        call output(' (2) nLevels = ', advance='no')
        call output(L2gp2%nLevels, advance='yes')
        ShapesDontMatch = .true.
      endif
      if ( L2gp1%nFreqs /= L2gp2%nFreqs ) then
        call output(' (1) nFreqs = ', advance='no')
        call output(L2gp1%nFreqs, advance='yes')
        call output(' (2) nFreqs = ', advance='no')
        call output(L2gp2%nFreqs, advance='yes')
        ShapesDontMatch = .true.
      endif
      if ( myDetails < 0 ) return
      if ( ShapesDontMatch ) then
        call output('Skipping further details because shapes dont match',&
          & advance='yes')
        return
      endif
      if ( any(l2gp1%pressures /= l2gp2%pressures)) then
          call dump ( l2gp1%pressures - l2gp2%pressures, &
            & 'l2gp%pressures (diff)', stats=stats, rms=rms )
      endif
      if ( any(l2gp1%latitude /= l2gp2%latitude)) then
          call dump ( l2gp1%latitude - l2gp2%latitude, &
            & 'l2gp%latitude (diff)', stats=stats, rms=rms )
      endif
      if ( any(l2gp1%longitude /= l2gp2%longitude)) then
          call dump ( l2gp1%longitude - l2gp2%longitude, &
            & 'l2gp%longitude (diff)', stats=stats, rms=rms )
      endif
      if ( any(l2gp1%solarTime /= l2gp2%solarTime)) then
        call dump ( l2gp1%solarTime - l2gp2%solarTime, &
          & 'l2gp%solarTime (diff)', stats=stats, rms=rms )
      endif
      if ( any(l2gp1%solarZenith /= l2gp2%solarZenith)) then
        call dump ( l2gp1%solarZenith - l2gp2%solarZenith, &
          & 'l2gp%solarZenith (diff)', stats=stats, rms=rms )
        badChunks = .true.
      endif
      if ( any(l2gp1%losAngle /= l2gp2%losAngle)) then
        call dump ( l2gp1%losAngle - l2gp2%losAngle, &
          & 'l2gp%losAngle (diff)', stats=stats, rms=rms )
      endif
      if ( any(l2gp1%geodAngle /= l2gp2%geodAngle)) then
        call dump ( l2gp1%geodAngle - l2gp2%geodAngle, &
          & 'l2gp%geodAngle (diff)', stats=stats, rms=rms )
        badChunks = .true.
      endif
      if ( any(l2gp1%time /= l2gp2%time)) then
        call dump ( l2gp1%time - l2gp2%time, &
          & 'l2gp%time (diff)', stats=stats, rms=rms )
        badChunks = .true.
      endif
      if ( any(l2gp1%chunkNumber /= l2gp2%chunkNumber)) then
        call dump ( l2gp1%chunkNumber - l2gp2%chunkNumber, &
          & 'l2gp%chunkNumber (diff)', stats=stats, rms=rms )
      endif
      
      if ( associated(l2gp1%frequency) .and.  associated(l2gp2%frequency)) then
        if ( any(l2gp1%frequency /= l2gp2%frequency)) then
          call dump ( l2gp1%frequency - l2gp2%frequency, &
            & 'l2gp%frequency (diff)', stats=stats, rms=rms )
        endif
      endif
      
      if ( present ( ignoreBadChunks ) ) then
        badChunks = ignoreBadChunks .and. badChunks
      else
        badChunks = .false.
      endif
      if ( myDetails < 1 ) return
      ! This is a very bad idea--to redefine arrays in a procedure
      ! whose ostensible purpose is merely to dump them
      ! It would be smarter to allocate a temp array and use
      ! that, remembering to deallocate it before returning
      if ( badChunks ) then
        do instance=1, L2gp1%nTimes
          if ( l2gp1%geodAngle(instance) /= l2gp2%geodAngle(instance) &
            & .or. &
            & l2gp1%solarZenith(instance) /= l2gp2%solarZenith(instance) &
            & .or. &
            & l2gp1%time(instance) /= l2gp2%time(instance) &
            & ) then
            l2gp2%l2gpValue(:,:,instance) = l2gp1%l2gpValue(:,:,instance)
            l2gp2%l2gpPrecision(:,:,instance) = l2gp1%l2gpPrecision(:,:,instance)
            l2gp2%status(instance) = l2gp1%status(instance)
            l2gp2%quality(instance) = l2gp1%quality(instance)
            badInstances = badInstances + 1
          endif
        enddo
        call output('Number of bad instances of l2gp2 reset to l2gp1 ', advance='no')
        call output(badInstances, advance='yes')
      endif
      if ( any(l2gp1%l2gpValue /= l2gp2%l2gpValue)) then
        call dump ( real(l2gp1%l2gpValue - l2gp2%l2gpValue, r8), &
          & 'l2gp%l2gpValue (diff)', stats=stats, rms=rms )
          !& 'l2gp%l2gpValue (diff)', FillValue=FillValue )
      endif
      if ( any(l2gp1%l2gpPrecision /= l2gp2%l2gpPrecision)) then
        call dump ( real(l2gp1%l2gpPrecision - l2gp2%l2gpPrecision, r8), &
          & 'l2gp%l2gpPrecision (diff)', stats=stats, rms=rms )
      endif
      
      if ( any(l2gp1%status /= l2gp2%status)) then
        call dump (l2gp1%status - l2gp2%status, &
          & 'l2gp%status (diff)', stats=stats, rms=rms )
      endif
      if ( any(l2gp1%quality /= l2gp2%quality)) then
        call dump ( l2gp1%quality - l2gp2%quality, &
          & 'l2gp%quality (diff)', stats=stats, rms=rms )
      endif
      
  end subroutine DiffL2GPData
    
  ! ------------------------------------------ DiffL2GPFiles ------------
  subroutine DiffL2GPFiles ( file1, file2, &
    & Details, wholeArray, stats, rms, ignoreBadChunks, swList, showMissing )
    ! Show diff between swaths in file1 and file2 down to level of Details
    ! Dummy arguments
    character (len=*), intent(in) :: file1 ! Name of file 1
    character (len=*), intent(in) :: file2 ! Name of file 2
    integer, intent(in), optional :: DETAILS ! <=0 => Don't diff data fields
    !                                        ! -1 Skip even geolocation fields
    !                                        ! -2 Skip all but name
    !                                        ! >0 Diff even data fields
    !                                        ! Default 1
    !
    ! The following parameters, if present, will override Details
    logical, intent(in), optional :: WHOLEARRAY
    logical, intent(in), optional :: RMS   ! if TRUE, just print mean, rms
    logical, intent(in), optional :: STATS   ! if TRUE, just print stats
    logical, intent(in), optional :: IGNOREBADCHUNKS
    ! swList currently used only if showMissing is TRUE
    character (len=*), optional, intent(in) :: swList
    logical, intent(in), optional :: showMissing   ! if TRUE, just show which
    ! Local                                         swaths are missing from other
    logical, parameter            :: countEmpty = .true.
    integer :: File1Handle
    integer :: File2Handle
    integer :: record_length
    integer :: i
    integer :: status
    integer :: the_hdfVersion1
    integer :: the_hdfVersion2
    logical :: file_exists
    integer :: file_access
    integer :: listsize
    integer :: noSwaths
    integer :: noSwaths2
    integer :: noUnique
    logical :: myShowMissing
    character (len=MAXSWATHNAMESBUFSIZE) :: swathList1
    character (len=MAXSWATHNAMESBUFSIZE) :: swathList2
    character (len=MAXSWATHNAMESBUFSIZE) :: swathUnique
    character (len=L2GPNameLen) :: swath
    type (L2GPData_T) :: l2gp1
    type (L2GPData_T) :: l2gp2
    ! Executable code
    myShowMissing = .false.
    if ( present(showMissing) ) myShowMissing=showMissing
    file_exists = ( mls_exists(trim(File1)) == 0 )
    if ( .not. file_exists ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'File 1 not found; make sure the name and path are correct' &
        & // trim(file1) )
    endif
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
    ! Are we merely to point out which swaths are missing from the other file?
    if ( myShowMissing ) then
      if ( present(swList) ) then
        call GetUniqueList(swList, swathUnique, noUnique, countEmpty, &
          & str2=trim(swathList1))
        if ( noUnique > 0 ) then
          call output('swaths missing from ' // trim(File1), advance='yes')
          call output(trim(swathUnique), advance='yes')
        else
          call output('All swaths in ' // trim(File1), advance='yes')
        endif
        call GetUniqueList(swList, swathUnique, noUnique, countEmpty, &
          & str2=trim(swathList2))
        if ( noUnique > 0 ) then
          call output('swaths missing from ' // trim(File2), advance='yes')
          call output(trim(swathUnique), advance='yes')
        else
          call output('All swaths in ' // trim(File2), advance='yes')
        endif
      else
        call output('Comparing swaths in ' // trim(File1), advance='no')
        call output(' with ' // trim(File2), advance='yes')
        call GetUniqueList(trim(swathList1), swathUnique, noUnique, countEmpty, &
          & str2=trim(swathList2))
         ! print *, 'swathList1: ', trim(swathList1)
         ! print *, 'swathList2: ', trim(swathList2)
         ! print *, 'list1 not in list2 : ', noUnique
        if ( noUnique > 0 ) then
          call output('swaths only in ' // trim(File1), advance='yes')
          call output(trim(swathUnique), advance='yes')
        endif
        call GetUniqueList(trim(swathList2), swathUnique, noUnique, countEmpty, &
          & str2=trim(swathList1))
        ! print *, 'list2 not in list1 : ', noUnique
        if ( noUnique > 0 ) then
          call output('swaths only in ' // trim(File2), advance='yes')
          call output(trim(swathUnique), advance='yes')
        endif
      endif
      return
    endif
    File1Handle = mls_io_gen_openF('swopen', .TRUE., status, &
       & record_length, DFACC_READ, FileName=File1, &
       & hdfVersion=the_hdfVersion1, debugOption=.false. )
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to open L2gp file: " // trim(File1) // ' for diff')
    File2Handle = mls_io_gen_openF('swopen', .TRUE., status, &
       & record_length, DFACC_READ, FileName=File2, &
       & hdfVersion=the_hdfVersion2, debugOption=.false. )
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to open L2gp file: " // trim(File2) // ' for diff')

    ! Loop over swaths in file 1

    ! print *, 'swathList1: ', trim(swathList1)
    do i = 1, noSwaths
      call GetStringElement (trim(swathList1), swath, i, countEmpty )
      if ( len_trim(swath) < 1 ) then
        call output('(Ignoring blank swath name in ' // &
          &  trim(File1), advance='yes')
        cycle
      endif
      status = stringElementNum(swathList2, trim(swath), countEmpty)
      if ( status < 1 ) then
        call output('Swath ' // trim(swath) // ' not found in ' // &
          & trim(File2), advance='yes')
        cycle
      endif
      call blanks( 22+len_trim(swath), fillChar='-', advance='yes')
      call output( '---- swath name: ' // trim(swath) // ' ----', advance='yes')
      call blanks( 22+len_trim(swath), fillChar='-', advance='yes')
      call ReadL2GPData ( File1Handle, trim(swath), l2gp1, &
           & hdfVersion=the_hdfVersion1 )
      call ReadL2GPData ( File2Handle, trim(swath), l2gp2, &
           & hdfVersion=the_hdfVersion2 )
      call DiffL2GPData(l2gp1, l2gp2, &
        & details=details, wholeArray=wholeArray, rms=rms, stats=stats, &
        & ignoreBadChunks=ignoreBadChunks)
      call DestroyL2GPContents ( l2gp1 )
      call DestroyL2GPContents ( l2gp2 )
    enddo
    status = mls_io_gen_closeF('swclose', File1Handle, FileName=File1, &
      & hdfVersion=the_hdfVersion1, debugOption=.false.)
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2gp file: " // trim(File1) // ' after diff')
    status = mls_io_gen_closeF('swclose', File2Handle, FileName=File2, &
      & hdfVersion=the_hdfVersion1, debugOption=.false.)
    if ( status /= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
       & "Unable to close L2gp file: " // trim(File2) // ' after diff')
  end subroutine DiffL2GPFiles
    
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

  ! ------------------------------------------ DUMP_L2GP_DATABASE ------------

  subroutine DUMP_L2GP_DATABASE ( L2gp, Name, ColumnsOnly, Details, Fields )

    ! Dummy arguments
    type (l2gpData_T), intent(in) ::          L2GP(:)
    character(len=*), intent(in), optional :: Name
    logical, intent(in), optional ::          ColumnsOnly ! if true, dump only with columns
    integer, intent(in), optional :: DETAILS
    character(len=*), intent(in), optional :: fields ! ,-separated list of names

    ! Local variables
    integer :: i
    call output ( '============ L2GP Data Base ============', advance='yes' )
    call output ( ' ', advance='yes' )
    if ( present(name) ) then
      call output ( 'L2GP Database name: ', advance='no' )
      call output ( name, advance='yes' )
    endif
    if ( size(l2gp) < 1 ) then
      call output ( '**** L2GP Database empty ****', advance='yes' )
      return
    endif
    do i = 1, size(l2gp)
      call dump(l2gp(i), ColumnsOnly, Details, Fields)
    end do
      
  end subroutine DUMP_L2GP_DATABASE

  ! ------------------------------------------ DUMP_L2GP ------------


  subroutine Dump_L2GP ( L2gp, ColumnsOnly, Details, Fields )
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

    ! Local variables
    integer :: ierr
    logical :: myColumnsOnly
    integer :: MYDETAILS
    real(r8) :: FillValue
    integer :: ChunkFillValue
    character(len=Len(DATA_FIELDS)+Len(GEO_FIELDS)) :: myFields
    logical :: show

    ! Executable code
    myDetails = 1
    if ( present(details) ) myDetails = details
    myFields = ' '
    if ( present(fields) ) then
      myFields = fields
      myDetails = 2
    endif
    
    if( present(ColumnsOnly)) then
      myColumnsOnly = ColumnsOnly
    else
      myColumnsOnly = .false.
    endif

    if ( myColumnsOnly .and. l2gp%nLevels > 1 ) return
    
    FillValue = real(l2gp%MissingValue, r8)
    ChunkFillValue = int(l2gp%MissingValue)
    if ( showMe(.true., myFields, 'swathname') ) then
      call output ( 'L2GP Data: (swath name) ')
      call output ( trim(l2gp%name) )
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
    if ( showMe(.true., myFields, 'quantitytype') ) then
      call output ( 'Quantity type: ')
      if ( l2gp%QuantityType > 0 ) then
        call display_string ( lit_indices(l2gp%QuantityType), &
          &             strip=.true., advance='yes' )
      else
        call output ( '(the QType Index was 0) ', advance='yes')
      endif
    endif
    !  if ( myDetails < -1 ) return
    if ( showMe(myDetails > -2, myFields, 'ntimes') ) then
      call output ( 'nTimes: ')
      call output ( l2gp%nTimes, 5)
      call output ( '  nLevels: ')
      call output ( l2gp%nLevels, 3)
      call output ( '  nFreqs: ')
      call output ( l2gp%nFreqs, 3, advance='yes')
      call output ( 'Fill/Missing Values: ')
      call output ( l2gp%MissingValue, advance='yes')
      call output ( 'Fill/Missing Status Field: ')
      call output ( l2gp%MissingStatus, advance='yes')
     endif
    
     ! if ( myDetails < 0 ) return
    if ( showMe(myDetails > -1, myFields, 'pressure') ) &
      & call dump ( l2gp%pressures, trim(l2gp%verticalCoordinate) // 's:' )
      
    if ( showMe(myDetails > -1, myFields, 'latitude') ) &
      & call dump ( l2gp%latitude, 'Latitude:' )
      
    if ( showMe(myDetails > -1, myFields, 'longitude') ) &
      & call dump ( l2gp%longitude, 'Longitude:' )
      
    if ( showMe(myDetails > -1, myFields, 'solartime') ) &
      & call dump ( l2gp%solarTime, 'SolarTime:' )
      
    if ( showMe(myDetails > -1, myFields, 'solarzenith') ) &
      & call dump ( l2gp%solarZenith, 'SolarZenith:' )
      
    if ( showMe(myDetails > -1, myFields, 'LOSAngle') ) &
      & call dump ( l2gp%losAngle, 'LOSAngle:' )
      
    if ( showMe(myDetails > -1, myFields, 'geodAngle') ) &
      & call dump ( l2gp%geodAngle, 'geodAngle:' )
      
    if ( showMe(myDetails > -1, myFields, 'time') ) &
      & call dump ( l2gp%time, 'Time:' )
      
    if ( showMe(myDetails > -1, myFields, 'chunkNumber') ) &
      & call dump ( l2gp%chunkNumber, 'ChunkNumber:' )
      
      if ( showMe(myDetails > -1, myFields, 'pressure') .and. &
        & associated(l2gp%frequency) ) &
        & call dump ( l2gp%frequency, 'Frequencies:' )
      
      ! if ( myDetails < 1 ) return
    if ( showMe(myDetails > 0, myFields, 'l2gpvalue') ) &
      & call dump ( real(l2gp%l2gpValue, r8), 'L2GPValue:', &
        & FillValue=FillValue )
      
    if ( showMe(myDetails > 0, myFields, 'l2gpprecision') ) &
      & call dump ( real(l2gp%l2gpPrecision, r8), 'L2GPPrecision:', &
        & FillValue=FillValue )
      
    if ( showMe(myDetails > 0, myFields, 'status') ) &
      & call dump ( l2gp%status, 'Status:' )
      
    if ( showMe(myDetails > 0, myFields, 'quality') ) &
      & call dump ( l2gp%quality, 'Quality:' )
      
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
        showMe = ( index(LowerCase(fields), LowerCase(trim(field))) > 0 )
      endif
    end function showMe
  end subroutine Dump_L2GP
    
  !----------------------------------------  DumpL2GP_attributes_hdf5  -----
  subroutine DumpL2GP_attributes_hdf5(l2FileHandle, l2gp, swathName)

  use HDFEOS5, only: HE5T_NATIVE_INT, HE5T_NATIVE_REAL, HE5T_NATIVE_DOUBLE, &
    & HE5_SWattach, HE5_SWdetach, MLS_charType
  use he5_swapi, only: he5_swrdattr, he5_swrdlattr
  use PCFHdr, only:  GlobalAttributes_T, he5_readglobalattr
    ! Brief description of subroutine
    ! This subroutine dumps the attributes for an l2gp
    ! These include
    ! Global attributes which are file level
    ! Swath level
    ! Geolocation field attributes
    ! Data field attributes
    ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp
    integer, intent(in) :: l2FileHandle ! From swopen
    character (len=*), optional, intent(IN) :: swathName ! Defaults->l2gp%name

    ! Variables
    type (GlobalAttributes_T)         :: gAttributes
    character(len=255)                :: ProcessLevel
    integer                           :: DayofYear
    double precision                  :: TAI93At0zOfGranule
    real(rgp), dimension(MAXNLEVELS)  :: pressures
    integer, dimension(1)             :: ibuf

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
    integer :: rgp_type
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
    
    ! Begin
    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    call output ( 'L2GP Attributes: (swath name) ')
    call output ( trim(name), advance='yes' )
    if ( rgp == r4 ) then
      rgp_type = HE5T_NATIVE_REAL
    elseif ( rgp == r8 ) then
      rgp_type = HE5T_NATIVE_DOUBLE
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & 'Attributes have unrecognized numeric data type; should be r4 or r8')
    endif
    call List2Array(GeolocationTitles, theTitles, .true.)
    call List2Array(GeolocationUnits, theUnits, .true.)

    ! - -   G l o b a l   A t t r i b u t e s   - -
    call output ( '(Global Attributes) ', advance='yes')
    call he5_readglobalattr(l2FileHandle, gAttributes, &
     & ProcessLevel, DayofYear, TAI93At0zOfGranule, status)
    if ( status /= 0 ) then
      call output ('No global attributes found in file', advance='yes')
    else
      call dump(gAttributes%orbNum, 'Orbit numbers')
      call dump(gAttributes%orbPeriod, 'Orbit Periods')
      call output ('InstrumentName: ' // trim( gAttributes%InstrumentName  ))
      call output ('Process level: ' // trim(  gAttributes%ProcessLevel    ))
      call output ('Input version: ' // trim(  gAttributes%inputVersion    ))
      call output ('PGE version: ' // trim(    gAttributes%PGEVersion      ))
      call output ('Start UTC: ' // trim(      gAttributes%StartUTC        ))
      call output ('End UTC: ' // trim(        gAttributes%EndUTC          ))
      call dump_int ( gAttributes%GranuleMonth, 'Granule month:' )
      call dump_int ( gAttributes%GranuleDay, 'Granule day:' )
      call dump_int ( gAttributes%GranuleYear, 'Granule year:' )
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
        call GetStringHashElement (GeolocationTitles, &
          & GeoUniqueFieldDefinition, trim(theTitles(field)), &
          & abbr_uniq_fdef, .false.)
        call GetStringHashElement (UniqueFieldDefKeys, &
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
    call GetStringHashElement (lowercase(Species), &
      & SpUniqueFieldDefinition, trim(lowercase(species_name)), &
      & abbr_uniq_fdef, .false.)
    call GetStringHashElement (UniqueFieldDefKeys, &
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

  ! ----------------------------------  GetGeolocUnits  -----
  function GetGeolocUnits ( Title ) result( dim_string )

  ! Given a dim name type, e.g. l_vmr,
  ! returns corresponding units as a character string

    ! Dummy arguments
    character(len=*), intent(in)        :: Title
    character(len=CHARATTRLEN)                   :: dim_string

    ! Executable code
    dim_string = 'none'
    select case (Title)                                       
    case ( 'Latitude' )  
      dim_string = 'degrees'
    case ( 'Longitude' )  
      dim_string = 'degrees'
    case ( 'SolarZenithAngle' )  
      dim_string = 'degrees'
    case ( 'Time' )  
      dim_string = 'seconds'
    case ( 'LocalSolarTime' )  
      dim_string = 'seconds'
    case ( 'LineOfSightAngle' )  
      dim_string = 'degrees'
    case ( 'OrbitGeodeticAngle' )  
      dim_string = 'degrees'
    case ( 'ChunkNumber' )
      ! No units for this dimension
    case ( 'Pressure' )  
      dim_string = 'hPa'
    case ( 'Frequency' )  
      dim_string = 'GHz'

    end select                                                       

  end function GetGeolocUnits


!=============================================================================
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

!=============================================================================
end module L2GPData
!=============================================================================

!
! $Log$
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
