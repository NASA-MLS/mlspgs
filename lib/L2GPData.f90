! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module L2GPData                 ! Creation, manipulation and I/O for L2GP Data
!=============================================================================
  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use DUMP_0, only: DUMP
  use Hdf, only: DFACC_READ, DFNT_CHAR8, DFNT_FLOAT32, DFNT_INT32, DFNT_FLOAT64
  use HDFEOS, only: SWATTACH, SWCREATE, SWDETACH, SWINQDIMS
  use Intrinsic ! "units" type literals, beginning with L_
  use MLSCommon, only: R4, R8
  use MLSFiles, only: FILENOTFOUND, HDFVERSION_4, HDFVERSION_5, &
    & MLS_HDF_VERSION, MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_Warning, MLSMSG_Debug
  use MLSStrings, only: Capitalize, ExtractSubString, GetStringHashElement, &
    & ints2Strings, list2array, lowercase, strings2Ints
  use OUTPUT_M, only: OUTPUT
  use PCFHdr, only: GA_VALUE_LENGTH
  use STRING_TABLE, only: DISPLAY_STRING

  implicit none

  private
  public :: L2GPData_T
  public :: L2GPNameLen
  public :: AddL2GPToDatabase,  DestroyL2GPContents,  DestroyL2GPDatabase, &
    & Dump, ExpandL2GPDataInPlace, AppendL2GPData, &
    & ReadL2GPData, SetupNewL2GPRecord,  WriteL2GPData

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       & "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character(len=*), parameter, private :: ModuleName = &
       & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  interface DUMP
    module procedure DUMP_L2GP
    module procedure DUMP_L2GP_DataBase
  end interface

  interface ReadL2GPData
    module procedure ReadL2GPData_fileID
    module procedure ReadL2GPData_fileName
  end interface

! interface my_swwrattr 
!   module procedure my_swwrattr_snglarray
!   module procedure my_swwrattr_char
! end interface

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
  integer, parameter :: L2GPDEFAULT_HDFVERSION = HDFVERSION_4

  ! r4 corresponds to sing. prec. :: same as stored in files
  integer, public, parameter :: rgp = r4

  integer, parameter :: CHARATTRLEN = GA_VALUE_LENGTH
  real, parameter    :: UNDEFINED_VALUE = -999.99 ! Same as %template%badvalue
  integer, parameter :: L2GPNameLen = 80
  integer, parameter :: NumGeolocFields = 10

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
  logical, parameter :: USEINTS4STRINGS = .false.  
  
  ! So far, the nameIndex component of the main data type is never set
  logical, parameter :: NAMEINDEXEVERSET = .false.  
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
     integer :: QUANTITYTYPE = 0         ! E.g., l_temperature

     ! Now the dimensions of the data

     integer :: nTimes          ! Total number of profiles
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

     ! Unfortunately, currently we are neither reading nor writing this field
     ! As explained in several places below, the cwrappers beneath have failed
     ! to observe the right number of memory places required to represent
     ! the fortran data.
     ! paw has an idea of how to do it eventually:
     ! it will require writing a separate module in which to read/write
     ! character-valued swath data fields, named something like 
     ! char_swdata_m.f90
     ! In this new module, he5_swrdfld and he5_swwrfld will simply
     ! be integer externals, preventing any type conversion
     ! before it gives their start addresses to the c-wrappers beneath
     ! Then just make sure that the right number of bytes get transferred
     character (len=1), pointer, dimension(:) :: status=>NULL()
     !                (status is a reserved word in F90)
     real (rgp), pointer, dimension(:) :: quality=>NULL()
     ! Both the above dimensioned (nTimes)

     ! The dimensions for the quantity (if, e.g., coming from l2cf)
     ! character(len=CHARATTRLEN)        :: DIM_Names = '' ! ','-separated
     ! character(len=CHARATTRLEN)        :: DIM_Units = '' ! ','-separated
     ! character(len=CHARATTRLEN)        :: VALUE_Units = '' 
     real (rgp)                        :: MissingValue = UNDEFINED_VALUE
  end type L2GPData_T

contains ! =====     Public Procedures     =============================

  !------------------------------------------  SetupNewL2GPRecord  -----
  subroutine SetupNewL2GPRecord ( l2gp, nFreqs, nLevels, nTimes, &
    & FillIn)

    ! This routine sets up the arrays for an l2gp datatype.

    ! Dummy arguments
    type (L2GPData_T), intent(inout)  :: l2gp
    integer, intent(in), optional :: nFreqs            ! Dimensions
    integer, intent(in), optional :: nLevels           ! Dimensions
    integer, intent(in), optional :: nTimes            ! Dimensions
    logical, intent(in), optional :: FillIn    ! Fill with MissingValue

    ! Local variables
    integer :: useNFreqs, useNLevels, useNTimes

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

    ! Store the dimensionality

    l2gp%nTimes = useNTimes
    l2gp%nLevels = useNLevels
    l2gp%nFreqs = useNFreqs

    ! But allocate to at least one for times, freqs
 
    useNTimes=MAX(useNTimes,1)
    useNLevels=MAX(useNLevels,1)
    useNFreqs=MAX(useNFreqs,1)    

    ! Allocate the vertical coordinate

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
    if ( .not. present(FillIn) ) return
    if ( .not. FillIn ) return
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
    l2gp%status = ' '
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
       firstProf, lastProf, hdfVersion, HMOT)
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
      call ReadL2GPData_hdf4(L2FileHandle, swathname, l2gp, my_hmot, &
        & numProfs, firstProf, lastProf)
    elseif (myhdfVersion /= HDFVERSION_5) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Unrecognized hdfVersion passed to ReadL2GPData" )
    else
      call ReadL2GPData_hdf5(L2FileHandle, swathname, l2gp, my_hmot, &
        & numProfs, firstProf, lastProf)
    endif
    !print*,"In readl2gpdata: first/last/read prof=",firstProf,&
    !  lastProf,numProfs
  end subroutine ReadL2GPData_fileID

  ! ---------------------- ReadL2GPData_fileName  -----------------------------

  subroutine ReadL2GPData_fileName(fileName, swathname, l2gp, numProfs, &
       firstProf, lastProf, hdfVersion, HMOT)
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

  ! ---------------------- ReadL2GPData_hdf4  -----------------------------

  subroutine ReadL2GPData_hdf4(L2FileHandle, swathname, l2gp, HMOT, &
    & numProfs, firstProf, lastProf)
  use MLSHDFEOS, only: mls_swdiminfo, mls_swrdfld
    !------------------------------------------------------------------------

    ! This routine reads an L2GP file, assuming it is hdfeos2 format
    ! (i.e. hdf4) returning a filled data structure and the
    ! number of profiles read.

    ! Arguments

    character (len=*), intent(in) :: swathname ! Name of swath
    integer, intent(in) :: L2FileHandle ! Returned by swopen
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( l2GPData_T ), intent(out) :: l2gp ! Result
    character, intent(in) :: HMOT   ! 'H' 'M'(def) 'O' 'T'
    integer, intent(out), optional :: numProfs ! Number actually read

    ! Local Parameters
    character (len=*), parameter :: SZ_ERR = 'Failed to get size of &
         &dimension '
    character (len=*), parameter :: MLSMSG_INPUT = 'Error in input argument '
    character (len=*), parameter :: MLSMSG_L2GPRead = 'Unable to read L2GP &
                                                     &field:'

    ! Local Variables
    character (len=80) :: DF_Name
    character (len=80) :: DF_Precision
    character (len=80) :: list
    character (len=480) :: msr
    integer :: alloc_err, first, freq, lev, nDims, size, swid, status
    integer :: start(3), stride(3), edge(3), dims(3)
    integer :: nFreqs, nLevels, nTimes, nFreqsOr1, nLevelsOr1, myNumProfs
    logical :: firstCheck, lastCheck

    real, allocatable :: realSurf(:), realProf(:), real3(:,:,:)
    logical :: dontfail
!  How to deal with status and columnTypes? swrfld fails
!  with char data on Linux
!  Have recourse to ints2Strings and strings2Ints if USEINTS4STRINGS
!    character (LEN=8), allocatable :: the_status_buffer(:)
!    character (LEN=L2GPNameLen), allocatable :: the_status_buffer(:)
    integer, allocatable, dimension(:,:) :: string_buffer

    ! Don't fail when trying to read an mls-specific field 
    ! if the file is from another Aura instrument
    dontfail = (HMOT /= 'M')

    ! Attach to the swath for reading

    l2gp%Name = swathname

    select case (HMOT)
    case ('H')
      swid = swattach(L2FileHandle, 'HIRDLS')
      DF_Name = TRIM(l2gp%Name)
      DF_Precision = TRIM(l2gp%Name) // 'Precision'
      l2gp%MissingValue = -999.
    case ('M')
      swid = swattach(L2FileHandle, TRIM(l2gp%Name))
      DF_Name = DATA_FIELD1
      DF_Precision = DATA_FIELD2
    case default
    end select

    if ( swid == -1) call MLSMessage ( MLSMSG_Error, ModuleName, 'Failed to &
         &attach to hdfeos2 swath interface for reading: ' // trim(swathname))

    ! Get dimension information

    lev = 0
    freq = 0

    nDims = swinqdims(swid, list, dims)
    if ( nDims == -1) call MLSMessage ( MLSMSG_Error, ModuleName, 'Failed to &
         &get dimension information on hdfeos2 swath ' // trim(swathname))
    if ( INDEX(list,'nLevels') /= 0 ) lev = 1

    select case (HMOT)
    case ('H')
      if ( INDEX(list,'Chan') /= 0 ) freq = 1
    case ('M')
      if ( INDEX(list,'Freq') /= 0 ) freq = 1
    case default
    end select

    size = mls_swdiminfo(swid, 'nTimes', hdfVersion=HDFVERSION_4)
    l2gp%nTimes = size
    nTimes=size
    if ( lev == 0 ) then
       nLevels = 0
    else
      size = mls_swdiminfo(swid, 'nLevels', hdfVersion=HDFVERSION_4)
      nLevels = size

    end if

    if ( freq == 1 .and. HMOT == 'M') then
      size = mls_swdiminfo(swid, 'nFreqs', hdfVersion=HDFVERSION_4)
      nFreqs = size
    elseif ( freq == 1 .and. HMOT == 'H') then
      size = mls_swdiminfo(swid, 'nChans', hdfVersion=HDFVERSION_4)
      nFreqs = size
    else
      nFreqs = 0
    end if

    ! Check optional input arguments

    firstCheck = present(firstProf)
    lastCheck = present(lastProf)

    if ( firstCheck ) then

       if ( (firstProf >= l2gp%nTimes) .OR. (firstProf < 0) ) then
          msr = MLSMSG_INPUT // 'firstProf ' // trim(swathname)
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       else
          first = firstProf
       end if

    else

       first = 0

    end if

    if ( lastCheck ) then

       if ( lastProf < first ) then
          msr = MLSMSG_INPUT // 'lastProf ' // trim(swathname)
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

       if ( lastProf >= nTimes ) then
          myNumProfs = nTimes - first
       else
          myNumProfs = lastProf - first + 1
       end if

    else

       myNumProfs = nTimes - first

    end if

    ! Allocate result

    call SetupNewL2GPRecord ( l2gp, nFreqs=nFreqs, nLevels=nLevels, &
    & nTimes=myNumProfs, FillIn = .true. )

    ! Allocate temporary arrays

    nFreqsOr1=MAX(nFreqs,1)
    nLevelsOr1=MAX(nLevels, 1)
    allocate ( realProf(myNumProfs), realSurf(l2gp%nLevels), &
      &   string_buffer(1,myNumProfs), &
      &   real3(nFreqsOr1,nLevelsOr1,myNumProfs), STAT=alloc_err )
    if ( alloc_err /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//' various things in ReadL2GPData' )

    ! Read the horizontal geolocation fields

    start(1) = 0
    start(2) = 0
    start(3) = first
    stride = 1
    edge(1) = nFreqsOr1
    edge(2) = nLevelsOr1
    edge(3) = myNumProfs
    status=0

    status = mls_swrdfld(swid, 'Latitude', start(3:3), stride(3:3), edge(3:3), &
      &    realProf, hdfVersion=HDFVERSION_4)
    l2gp%latitude = realProf

    status = mls_swrdfld(swid, 'Longitude', start(3:3), stride(3:3), edge(3:3), &
      &    realProf, hdfVersion=HDFVERSION_4)
    l2gp%longitude = realProf

    status = mls_swrdfld(swid, 'Time', start(3:3), stride(3:3), edge(3:3), &
      &    l2gp%time, hdfVersion=HDFVERSION_4)

    status = mls_swrdfld(swid, 'LocalSolarTime', start(3:3), stride(3:3), edge(3:3), &
      &    realProf, hdfVersion=HDFVERSION_4)
    l2gp%solarTime = realProf

    status = mls_swrdfld(swid, 'SolarZenithAngle', start(3:3), stride(3:3), edge(3:3), &
      &    realProf, hdfVersion=HDFVERSION_4)
    l2gp%solarZenith = realProf

    status = mls_swrdfld(swid, 'LineOfSightAngle', start(3:3), stride(3:3), edge(3:3), &
      &    realProf, hdfVersion=HDFVERSION_4, dontfail=dontfail)
    l2gp%losAngle = realProf

    status = mls_swrdfld(swid, 'OrbitGeodeticAngle', start(3:3), stride(3:3), edge(3:3), &
      &    realProf, hdfVersion=HDFVERSION_4, dontfail=dontfail)
    l2gp%geodAngle = realProf

    status = mls_swrdfld(swid, 'ChunkNumber', start(3:3), stride(3:3), edge(3:3), &
      &    l2gp%chunkNumber, hdfVersion=HDFVERSION_4, dontfail=dontfail)

    ! Read the pressures vertical geolocation field, if it exists

    if ( lev /= 0 ) then

       status = mls_swrdfld(swid, 'Pressure', start(2:2), stride(2:2), edge(2:2), &
         & realSurf, hdfVersion=HDFVERSION_4, dontfail=dontfail)
       l2gp%pressures = realSurf

    end if

    ! Read the frequency geolocation field, if it exists

    if ( freq == 1 ) then

       edge(1) = l2gp%nFreqs

       status = mls_swrdfld(swid, 'Frequency', start(1:1), stride(1:1), edge(1:1), &
         & l2gp%frequency, hdfVersion=HDFVERSION_4)

    end if

    ! Read the data fields that may have 1-3 dimensions

    if ( freq == 1 ) then

       status = mls_swrdfld(swid, trim(DF_Name), start, stride, edge, real3, &
         & hdfVersion=HDFVERSION_4)
       l2gp%l2gpValue = real3

       status = mls_swrdfld(swid, trim(DF_Precision), start, stride, edge, real3, &
         & hdfVersion=HDFVERSION_4)
       l2gp%l2gpPrecision = real3

    else if ( lev == 1 ) then

      status = mls_swrdfld( swid, trim(DF_Name), start(2:3), stride(2:3), &
        &   edge(2:3), real3(1,:,:), hdfVersion=HDFVERSION_4 )
      l2gp%l2gpValue = real3
      
      status = mls_swrdfld( swid, trim(DF_Precision), start(2:3), stride(2:3), &
        & edge(2:3), real3(1,:,:), hdfVersion=HDFVERSION_4 )
      l2gp%l2gpPrecision = real3
      
    else

       status = mls_swrdfld( swid, trim(DF_Name), start(3:3), stride(3:3), edge(3:3), &
         &   real3(1,1,:), hdfVersion=HDFVERSION_4 )
       l2gp%l2gpValue = real3

       status = mls_swrdfld( swid, trim(DF_Precision), start(3:3), stride(3:3), edge(3:3), &
            real3(1,1,:), hdfVersion=HDFVERSION_4 )
       l2gp%l2gpPrecision = real3

    end if

    ! Read the data fields that are 1-dimensional

!??? There appears to be a problem with reading character data using
!??? HDF-EOS.  HDF_EOS's swrdfld expects a C void* pointer, no matter what
!??? the type of object.  Burkhard Burow's cfortran macros need to be told
!??? whether the argument is character, because compilers represent character
!??? arguments in various ways, usually different from non-character arguments.
!??? I.e., never the twain shall meet.
!    status = swrdfld(swid, DATA_FIELD3,start(3:3),stride(3:3),edge(3:3),&
!      l2gp%status)
!    status = swrdfld(swid, DATA_FIELD3,start(3:3),stride(3:3),edge(3:3),&
!      the_status_buffer)
! These lines commented out as they make NAG core dump on the deallocate statement.
! below.
!    if ( status == -1 ) then
!      msr = MLSMSG_L2GPRead // DATA_FIELD3
!      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
!    end if
!    l2gp%status = the_status_buffer(:)(1:1)

    ! (   see note above concerning char_swdata_m.f90   )

    l2gp%status = ' ' ! So it has a value.

!   if(USEINTS4STRINGS) then
!       status = mls_swrdfld(swid, DATA_FIELD3,start(3:3),stride(3:3),edge(3:3),&
!         string_buffer, hdfVersion=HDFVERSION_4)
!       call ints2Strings(string_buffer, l2gp%status)
!    end if

    status = mls_swrdfld(swid, 'Quality', start(3:3), stride(3:3), edge(3:3), &
      &   realProf, hdfVersion=HDFVERSION_4, dontfail=dontfail)
    l2gp%quality = realProf

    ! Deallocate local variables

    deallocate ( realSurf, realProf, real3, STAT=alloc_err )
    if ( alloc_err /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         'Failed deallocation of local real variables.' )

    deallocate ( string_buffer, STAT=alloc_err )
    if ( alloc_err /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         'Failed deallocation of status buffer.' )

    !  After reading, detach from swath interface

    status = swdetach(swid)
    if ( status == -1) call MLSMessage ( MLSMSG_Error, ModuleName, &
         & 'Failed to detach from swath interface after reading.' )

    ! Set numProfs if wanted
    if ( present(numProfs) ) numProfs=myNumProfs

    !-----------------------------
  end subroutine ReadL2GPData_hdf4
  !-----------------------------

  ! ------------------- ReadL2GPData_hdf5 ----------------

  subroutine ReadL2GPData_hdf5(L2FileHandle, swathname, l2gp, HMOT, &
    & numProfs, firstProf, lastProf)
  use HDFEOS5, only: HE5_swattach, HE5_swdetach, HE5_SWINQDIMS
  use MLSHDFEOS, only: mls_swdiminfo, mls_swrdfld
    !------------------------------------------------------------------------

    ! This routine reads an L2GP file, returning a filled data structure and the !
    ! number of profiles read.

    ! Arguments

    character (LEN=*), intent(IN) :: swathname ! Name of swath
    integer, intent(IN) :: L2FileHandle ! Returned by HE5_swopen
    integer, intent(IN), optional :: firstProf, lastProf ! Defaults to first and last
    type( L2GPData_T ), intent(OUT) :: l2gp ! Result
    character, intent(in) :: HMOT   ! 'H' 'M'(def) 'O' 'T'
    integer, intent(OUT),optional :: numProfs ! Number actually read

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
    character (LEN=480) :: msr

    integer :: alloc_err, first, freq, lev, nDims, size, swid, status
    integer :: start(3), stride(3), edge(3), dims(3)
    integer :: nFreqs, nLevels, nTimes, nFreqsOr1, nLevelsOr1, myNumProfs
    logical :: firstCheck, lastCheck, timeIsUnlim

    real, allocatable :: realFreq(:), realSurf(:), realProf(:), real3(:,:,:)
    logical :: dontfail
!  How to deal with status and columnTypes? swrfld fails
!  with char data on Linux with HDF4. With HDF5 we may or may not need to
!  Have recourse to ints2Strings and strings2Ints if USEINTS4STRINGS
!    character (LEN=8), allocatable :: the_status_buffer(:)
!    character (LEN=L2GPNameLen), allocatable :: the_status_buffer(:)
    integer, allocatable, dimension(:,:) :: string_buffer
    
    ! Don't fail when trying to read an mls-specific field 
    ! if the file is from another Aura instrument
    dontfail = (HMOT /= 'M')

    ! Attach to the swath for reading
    !print*," in readl2gpdata_hdf5: first/last=",firstprof,lastprof
    l2gp%Name = swathname
    
    select case (HMOT)
    case ('H')
      swid = HE5_SWattach(L2FileHandle, 'HIRDLS')
      DF_Name = TRIM(l2gp%Name)
      DF_Precision = TRIM(l2gp%Name) // 'Precision'
      l2gp%MissingValue = -999.
    case ('M')
      swid = HE5_SWattach(L2FileHandle, l2gp%Name)
      DF_Name = DATA_FIELD1
      DF_Precision = DATA_FIELD2
    case default
    end select
    if (swid == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &attach to hdfeos5 swath interface for reading' // trim(swathname))

    ! Get dimension information

    lev = 0
    freq = 0

    nDims = HE5_SWinqdims(swid, list, dims)
    !print*,"just called inqdims: nDims=",ndims,"list=",list,"dims=",dims
    if (nDims == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &get dimension information on hdfeos5 swath ' // trim(swathname))
    if ( index(list,'nLevels') /= 0 ) lev = 1
    if ( index(list,'Freq') /= 0 ) freq = 1
!    if ( index(list,'Unlim') /= 0 ) then 
!      timeIsUnlim = .TRUE.
!    else
!      timeIsUnlim = .FALSE.
!    endif

    size = mls_swdiminfo(swid, 'nTimes', hdfVersion=HDFVERSION_5)
    ! This will be wrong if timeIsUnlim .eq. .TRUE. . 
    ! HE5_SWdiminfo returns 1 instead of the right answer.
      
    l2gp%nTimes = size
    nTimes=size

    if (lev == 0) then
       nLevels = 0
    else
      size = mls_swdiminfo(swid, 'nLevels', hdfVersion=HDFVERSION_5)
       nLevels = size

    endif

    if ( freq == 1 .and. HMOT == 'M') then
      size = mls_swdiminfo(swid, 'nFreqs', hdfVersion=HDFVERSION_5)
       nFreqs = size
    elseif ( freq == 1 .and. HMOT == 'H') then
      size = mls_swdiminfo(swid, 'nChans', hdfVersion=HDFVERSION_5)
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

       myNumProfs = nTimes - first

    endif

    ! Allocate result

    call SetupNewL2GPRecord (l2gp, nFreqs=nFreqs, nLevels=nLevels, &
      &  nTimes=mynumProfs, FillIn = .true.)

    ! Allocate temporary arrays

    nFreqsOr1=max(nFreqs,1)
    nLevelsOr1=max(nLevels, 1)
    allocate(realProf(myNumProfs), realSurf(l2gp%nLevels), &
      &   realFreq(l2gp%nFreqs), &
      &   string_buffer(1,myNumProfs), &
      &   real3(nFreqsOr1,nLevelsOr1,myNumProfs), STAT=alloc_err)
    if ( alloc_err /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//' various things in ReadL2GPData' )

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
      edge(3:3), realProf, hdfVersion=HDFVERSION_5)
    l2gp%latitude = realProf

    status = mls_SWrdfld(swid, 'Longitude', start(3:3), stride(3:3), edge(3:3),&
      &    realProf, hdfVersion=HDFVERSION_5)
    l2gp%longitude = realProf

    status = mls_SWrdfld(swid, 'Time', start(3:3), stride(3:3), edge(3:3),&
      &    l2gp%time, hdfVersion=HDFVERSION_5)

    status = mls_SWrdfld(swid, 'LocalSolarTime', start(3:3), stride(3:3), edge(3:3),&
      &    realProf, hdfVersion=HDFVERSION_5)
    l2gp%solarTime = realProf

    status = mls_SWrdfld(swid, 'SolarZenithAngle', start(3:3), stride(3:3), edge(3:3),&
      &    realProf, hdfVersion=HDFVERSION_5)
    l2gp%solarZenith = realProf

    status = mls_SWrdfld(swid, 'LineOfSightAngle', start(3:3), stride(3:3), edge(3:3),&
      &    realProf, hdfVersion=HDFVERSION_5, dontfail=dontfail)
    l2gp%losAngle = realProf

    status = mls_SWrdfld(swid, 'OrbitGeodeticAngle', start(3:3), stride(3:3), edge(3:3),&
      &   realProf, hdfVersion=HDFVERSION_5, dontfail=dontfail)
    l2gp%geodAngle = realProf

    status = mls_SWrdfld(swid, 'ChunkNumber', start(3:3), stride(3:3), edge(3:3),&
      &    l2gp%chunkNumber, hdfVersion=HDFVERSION_5, dontfail=dontfail)

    ! Read the pressures vertical geolocation field, if it exists

    if (lev /= 0) then

       status = mls_SWrdfld(swid,'Pressure',start(2:2),stride(2:2), edge(2:2),&
         & realSurf, hdfVersion=HDFVERSION_5, dontfail=dontfail)
       l2gp%pressures = realSurf

    endif

    ! Read the frequency geolocation field, if it exists

    if (freq == 1) then

       edge(1) = l2gp%nFreqs

       status = mls_SWrdfld(swid,'Frequency',start(1:1),stride(1:1),edge(1:1),&
         & realFreq, hdfVersion=HDFVERSION_5, dontfail=dontfail)
       l2gp%frequency = realFreq

    endif

    ! Read the data fields that may have 1-3 dimensions

    if ( freq == 1) then

       status = mls_SWrdfld(swid, trim(DF_Name), start, stride, edge, real3, &
         & hdfVersion=HDFVERSION_5)
       l2gp%l2gpValue = real3

       status = mls_SWrdfld(swid, trim(DF_Precision), start, stride, edge, real3, &
         & hdfVersion=HDFVERSION_5)
       l2gp%l2gpPrecision = real3

    else if ( lev == 1) then

       status = mls_SWrdfld( swid, trim(DF_Name), start(2:3), stride(2:3), &
            edge(2:3), real3(1,:,:), hdfVersion=HDFVERSION_5 )
       l2gp%l2gpValue = real3

       status = mls_SWrdfld( swid, trim(DF_Precision), start(2:3), stride(2:3), &
            edge(2:3), real3(1,:,:), hdfVersion=HDFVERSION_5 )
       l2gp%l2gpPrecision = real3

    else

       status = mls_SWrdfld(swid,trim(DF_Name),start(3:3),stride(3:3),edge(3:3),&
         &   real3(1,1,:), hdfVersion=HDFVERSION_5 )
       l2gp%l2gpValue = real3

       status = mls_SWrdfld( swid, trim(DF_Precision), start(3:3), stride(3:3), edge(3:3), &
            real3(1,1,:), hdfVersion=HDFVERSION_5 )
       l2gp%l2gpPrecision = real3

    endif
    
    ! Read the data fields that are 1-dimensional

!??? There appears to be a problem with reading character data using
!??? HDF-EOS.  HDF_EOS's swrdfld expects a C void* pointer, no matter what
!??? the type of object.  Burkhard Burow's cfortran macros need to be told
!??? whether the argument is character, because compilers represent character
!??? arguments in various ways, usually different from non-character arguments.
!??? I.e., never the twain shall meet.
!    status = swrdfld(swid, DATA_FIELD3,start(3:3),stride(3:3),edge(3:3),&
!      l2gp%status)
!    status = swrdfld(swid, DATA_FIELD3,start(3:3),stride(3:3),edge(3:3),&
!      the_status_buffer)
! These lines commented out as they make NAG core dump on the deallocate statement.
! below.
!    if ( status == -1 ) then
!      msr = MLSMSG_L2GPRead // DATA_FIELD3
!      call MLSMessage ( MLSMSG_Error, ModuleName, msr )
!    end if
!    l2gp%status = the_status_buffer(:)(1:1)
    
    !   The above note was copied direct from the HDF4 version. The HDF5
    ! version has similar problems so these lines are commented too.
    !         status = HE5_SWrdfld(swid, DATA_FIELD3,start(3:3),&
    !    stride(3:3),edge(3:3), l2gp%status)
    !
    ! (   see note above concerning char_swdata_m.f90   )

    l2gp%status = ' ' ! So it has a value.

   if(USEINTS4STRINGS) then
!       status = mls_swrdfld(swid,DATA_FIELD3,start(3:3),stride(3:3),edge(3:3),&
!         string_buffer, hdfVersion=HDFVERSION_5)
!       call ints2Strings(string_buffer, l2gp%status)
    else
      call MLSMessage(MLSMSG_Debug, ModuleName, &
        "reading of status field disabled")
      status=0
    end if

    status = mls_SWrdfld(swid, 'Quality', start(3:3), stride(3:3),&
      edge(3:3),realProf, hdfVersion=HDFVERSION_5, dontfail=dontfail)
    l2gp%quality = realProf

    ! Deallocate local variables

    deallocate(realFreq, realSurf, realProf, real3, STAT=alloc_err)
    if ( alloc_err /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
         'Failed deallocation of local real variables.')

    deallocate ( string_buffer, STAT=alloc_err )
    if ( alloc_err /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         'Failed deallocation of status buffer.' )

    !  After reading, detach from HE5_SWath interface

    status = HE5_SWdetach(swid)
    if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &detach from swath interface after reading.')
    !print*," leaving readl2gpdata_hdf5: first/last/read=",&
    !  firstprof,lastprof,myNumProfs
    ! Set numProfs if wanted
    if (present(numProfs)) numProfs=myNumProfs


    !-----------------------------
  end subroutine ReadL2GPData_hdf5
  !-----------------------------

  ! --------------------------------------  OutputL2GP_createFile_hdf4  -----
  subroutine OutputL2GP_createFile_hdf4 (l2gp, L2FileHandle, swathName)
  use MLSHDFEOS, only: mls_dfldsetup, mls_gfldsetup, mls_swdefdim

    ! Brief description of subroutine
    ! This subroutine sets up the structural definitions in an empty L2GP file.

    ! Arguments

    integer, intent(in) :: L2FileHandle ! From swopen
    type( l2GPData_T ), intent(inout) :: l2gp
    character (len=*), optional, intent(in) :: swathName ! Defaults to l2gp%swathName

    ! Parameters

    character (len=*), parameter :: DIM_ERR = 'Failed to define dimension '
    character (len=*), parameter :: GEO_ERR = &
         & 'failed to define geolocation field '
    character (len=*), parameter :: DAT_ERR = 'Failed to define data field '

    ! Variables

    character (len=480) :: MSR
    character (len=132) :: NAME   ! From l2gp%name

    integer :: SWID, STATUS

    if ( present(swathName) ) then
       name=swathName
    else
       name=l2gp%name
    end if

    ! Create the swath within the file

    swid = swcreate(L2FileHandle, TRIM(name))
    if ( swid == -1 ) then
       msr = 'Failed to create swath ' // TRIM(name) &
        & // ' (maybe has the same name as another swath in this file?)'
    end if

    ! Define dimensions

    status = mls_swdefdim(swid, 'nTimes', l2gp%nTimes, &
      & hdfVersion=HDFVERSION_4)

    if ( l2gp%nLevels > 0 ) then
      status = mls_swdefdim(swid, 'nLevels', l2gp%nLevels, &
        & hdfVersion=HDFVERSION_4)
    end if

    if ( l2gp%nFreqs > 0 ) then
      status = mls_swdefdim(swid, 'nFreqs', l2gp%nFreqs, &
        & hdfVersion=HDFVERSION_4)
    end if

    ! Define horizontal geolocation fields using above dimensions

    status = mls_gfldsetup(swid, 'Latitude', 'nTimes', DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)
         
    status = mls_gfldsetup(swid, 'Longitude', 'nTimes', DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)

    status = mls_gfldsetup(swid, 'Time', 'nTimes', DEFAULTMAXDIM, &
      & DFNT_FLOAT64, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)

    status = mls_gfldsetup(swid, 'LocalSolarTime', 'nTimes', DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)

    status = mls_gfldsetup(swid, 'SolarZenithAngle', 'nTimes', DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)

    status = mls_gfldsetup(swid, 'LineOfSightAngle', 'nTimes', DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)

    status = mls_gfldsetup(swid, 'OrbitGeodeticAngle', 'nTimes', DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)

    status = mls_gfldsetup(swid, 'ChunkNumber', 'nTimes', DEFAULTMAXDIM, &
      & DFNT_INT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)

    if ( l2gp%nLevels > 0 ) then
      status = mls_gfldsetup(swid, 'Pressure', 'nLevels', DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)
    end if

    if ( l2gp%nFreqs > 0 ) then
      status = mls_gfldsetup(swid, 'Frequency', 'nFreqs', DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)
    end if

    ! Define data fields using above dimensions

    if ( (l2gp%nFreqs > 0) .AND. (l2gp%nLevels > 0) ) then

      status = mls_dfldsetup(swid, 'L2gpValue', 'nFreqs,nLevels,nTimes', &
      & DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)

      status = mls_dfldsetup(swid, 'L2gpPrecision', 'nFreqs,nLevels,nTimes', &
      & DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)

    else if ( l2gp%nLevels > 0 ) then

      status = mls_dfldsetup(swid, 'L2gpValue', 'nLevels,nTimes', &
      & DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)

      status = mls_dfldsetup(swid, 'L2gpPrecision', 'nLevels,nTimes', &
      & DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)

    else

      status = mls_dfldsetup(swid, 'L2gpValue', 'nTimes', &
      & DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)

      status = mls_dfldsetup(swid, 'L2gpPrecision', 'nTimes', &
      & DEFAULTMAXDIM, &
      & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
      & hdfVersion=HDFVERSION_4)

    end if

    status = mls_dfldsetup(swid, 'Status', 'nTimes', &
    & DEFAULTMAXDIM, &
    & DFNT_CHAR8, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
    & hdfVersion=HDFVERSION_4)

    status = mls_dfldsetup(swid, 'Quality', 'nTimes', &
    & DEFAULTMAXDIM, &
    & DFNT_FLOAT32, HDFE_NOMERGE, DEFAULT_CHUNKRANK, DEFAULT_CHUNKDIMS, &
    & hdfVersion=HDFVERSION_4)

    ! Detach from the swath interface.  This stores the swath info within the
    ! file and must be done before writing or reading data to or from the
    ! swath.

    status = swdetach(swid)
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Failed to detach from swath interface after definition.' )
    end if

    !--------------------------------------
  end subroutine OutputL2GP_createFile_hdf4
  !--------------------------------------

  !-----------------------------------------  OutputL2GP_writeGeo_hdf4  -----
  subroutine OutputL2GP_writeGeo_hdf4 (l2gp, l2FileHandle, swathName, offset)
  use MLSHDFEOS, only: mls_swwrfld

    ! Brief description of subroutine
    ! This subroutine writes the geolocation fields to an L2GP output file.

    ! Arguments

    type( l2GPData_T ), intent(inout) :: l2gp
    integer, intent(in) :: l2FileHandle ! From swopen
    character (len=*), intent(in), optional :: swathName ! Defaults to l2gp%name
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

    if ( present(swathName) ) then
       name=swathName
    else
       name=l2gp%name
    end if

    swid = swattach (l2FileHandle, name)

    ! Write data to the fields

    stride(1) = 1
    start(1) = 0     ! myOffset
    edge(1) = l2gp%nTimes

    status = mls_swwrfld(swid, 'Latitude', start, stride, edge, &
         real(l2gp%latitude), hdfVersion=HDFVERSION_4)

    status = mls_swwrfld(swid, 'Longitude', start, stride, edge, &
         real(l2gp%longitude), hdfVersion=HDFVERSION_4)

    status = mls_swwrfld(swid, 'Time', start, stride, edge, &
         l2gp%time, hdfVersion=HDFVERSION_4)

    status = mls_swwrfld(swid, 'LocalSolarTime', start, stride, edge, &
        real(l2gp%solarTime), hdfVersion=HDFVERSION_4)

    status = mls_swwrfld(swid, 'SolarZenithAngle', start, stride, edge, &
         real(l2gp%solarZenith), hdfVersion=HDFVERSION_4)

    status = mls_swwrfld(swid, 'LineOfSightAngle', start, stride, edge, &
         real(l2gp%losAngle), hdfVersion=HDFVERSION_4)

    status = mls_swwrfld(swid, 'OrbitGeodeticAngle', start, stride, edge, &
         real(l2gp%geodAngle), hdfVersion=HDFVERSION_4)

    status = mls_swwrfld(swid, 'ChunkNumber', start, stride, edge, &
         l2gp%chunkNumber, hdfVersion=HDFVERSION_4)

    if ( l2gp%nLevels > 0 ) then
       edge(1) = l2gp%nLevels
       ! start(1)=0 ! needed because offset may have made this /=0
       status = mls_swwrfld(swid, 'Pressure', start, stride, edge, &
            real(l2gp%pressures), hdfVersion=HDFVERSION_4)
    end if

    if ( l2gp%nFreqs > 0 ) then
       edge(1) = l2gp%nFreqs
       ! start(1)=0 ! needed because offset may have made this /=0
       status = mls_swwrfld(swid, 'Frequency', start, stride, edge, &
            real(l2gp%frequency), hdfVersion=HDFVERSION_4)
    end if

    ! Detach from the swath interface.  

    status = swdetach(swid)

    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach from swath interface' )
    end if

    !------------------------------------
  end subroutine OutputL2GP_writeGeo_hdf4
  !------------------------------------

  !----------------------------------------  OutputL2GP_writeData_hdf4  -----
  subroutine OutputL2GP_writeData_hdf4(l2gp, l2FileHandle, swathName, offset)
  use MLSHDFEOS, only: mls_swwrfld

    ! Brief description of subroutine
    ! This subroutine writes the data fields to an L2GP output file.

    ! Arguments

    type( l2GPData_T ), intent(inout) :: l2gp
    integer, intent(in) :: l2FileHandle ! From swopen
    character (len=*), intent(in), optional :: swathName ! Defaults to l2gp%name
    integer,intent(IN),optional::offset

    ! Parameters

    character (len=*), parameter :: WR_ERR = 'Failed to write data field '

    ! Variables

    character (len=132) :: name     ! Either swathName or l2gp%name

    integer :: status, myOffset
    integer :: start(3), stride(3), edge(3)
    integer :: swid

!  How to deal with status and columnTypes? swrfld fails
!  with char data on Linux
!  Have recourse to ints2Strings and strings2Ints if USEINTS4STRINGS
!    character (LEN=8), allocatable :: the_status_buffer(:)
!    character (LEN=L2GPNameLen), allocatable :: the_status_buffer(:)
    integer, allocatable, dimension(:,:) :: string_buffer

    ! Begin
    if (present(offset)) then
       myOffset=offset
    else
       myOffset=0
    endif
    if ( present(swathName) ) then
       name=swathName
    else
       name=l2gp%name
    end if
    ! Write data to the fields

    start = 0
    stride = 1
    ! start(3)= myOffset
    edge(1) = l2gp%nFreqs
    edge(2) = l2gp%nLevels
    edge(3) = l2gp%nTimes
    swid = swattach (l2FileHandle, name)
    if ( l2gp%nFreqs > 0 ) then
       ! Value and Precision are 3-D fields

       status = mls_swwrfld(swid, 'L2gpValue', start, stride, edge, &
            & RESHAPE(real(l2gp%l2gpValue), (/SIZE(l2gp%l2gpValue)/)), &
            & hdfVersion=HDFVERSION_4 )
       status = mls_swwrfld(swid, 'L2gpPrecision', start, stride, edge, &
            & RESHAPE(real(l2gp%l2gpPrecision), (/SIZE(l2gp%l2gpPrecision)/)), &
            & hdfVersion=HDFVERSION_4 )

    else if ( l2gp%nLevels > 0 ) then
       ! Value and Precision are 2-D fields

       status = mls_swwrfld( swid, 'L2gpValue', start(2:3), stride(2:3), &
            edge(2:3), real(l2gp%l2gpValue(1,:,:)), hdfVersion=HDFVERSION_4 )

       status = mls_swwrfld( swid, 'L2gpPrecision', start(2:3), stride(2:3), &
            edge(2:3), real(l2gp%l2gpPrecision(1,:,:) ), hdfVersion=HDFVERSION_4)
    else

       ! Value and Precision are 1-D fields
       status = mls_swwrfld( swid, 'L2gpValue', start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpValue(1,1,:) ), hdfVersion=HDFVERSION_4)
       status = mls_swwrfld( swid, 'L2gpPrecision', start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpPrecision(1,1,:) ), hdfVersion=HDFVERSION_4)
    end if

    ! 1-D status & quality fields

   if(USEINTS4STRINGS) then
!      allocate(string_buffer(1,l2gp%nTimes))
!      call strings2Ints(l2gp%status, string_buffer)
!      status = mls_swwrfld(swid, DATA_FIELD3, start(3:3), stride(3:3), edge(3:3), &
!           string_buffer, hdfVersion=HDFVERSION_4)
!      deallocate(string_buffer)
   else
!      status = mls_swwrfld(swid, DATA_FIELD3, start(3:3), stride(3:3), edge(3:3), &
!           l2gp%status, hdfVersion=HDFVERSION_4)
   end if
    !  l2gp%quality = 0 !??????? Why was this here !??? NJL
    status = mls_swwrfld(swid, 'Quality', start(3:3), stride(3:3), edge(3:3), &
         real(l2gp%quality), hdfVersion=HDFVERSION_4)

    !     Detach from the swath interface.

    status = swdetach(swid)
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach  from swath interface' )
    end if

    !-------------------------------------
  end subroutine OutputL2GP_writeData_hdf4
  !-------------------------------------

  ! --------------------------------------  OutputL2GP_createFile_hdf5  -----
  subroutine OutputL2GP_createFile_hdf5 (l2gp, L2FileHandle, swathName,nLevels)

  use HDFEOS5, only: HE5_SWcreate, HE5_SWdetach, &
    & HE5S_UNLIMITED_F, &
    & HE5T_NATIVE_CHAR, HE5T_NATIVE_DOUBLE, HE5T_NATIVE_INT, HE5T_NATIVE_FLOAT
  use MLSHDFEOS, ONLY : mls_dfldsetup, mls_gfldsetup, mls_swdefdim, &
    & MLS_SWSETFILL
    ! Brief description of subroutine
    ! This subroutine sets up the structural definitions in an empty L2GP file.

    ! Arguments

    integer, intent(in) :: L2FileHandle ! From swopen
    type( L2GPData_T ), intent(inout) :: l2gp
    character (LEN=*), optional, intent(IN) :: swathName ! Defaults to l2gp%swathName
    integer,optional::nLevels
    ! Parameters

    character (len=*), parameter :: DIM_ERR = 'Failed to define dimension '
    character (len=*), parameter :: GEO_ERR = &
         & 'Failed to define geolocation field '
    character (len=*), parameter :: DAT_ERR = 'Failed to define data field '

    ! Variables

    character (len=480) :: MSR
    character (len=132) :: NAME   ! From l2gp%name

    ! THESE ARE HDF5 CHUNKS, _NOT_ MLS ALONG-TRACK PROCESSING CHUNKS 
    integer,dimension(7)::CHUNK_DIMS
    integer::CHUNK_RANK
    integer::CHUNKTIMES,CHUNKFREQS,CHUNKLEVELS

    integer :: SWID, STATUS
    ! integer, external :: he5_swsetfill

    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    chunktimes=120      ! was 1
    chunkfreqs=1 ! better as nFreqs, but I have yet to see a case with nfreqs>1
    if(present(nLevels))then
       chunklevels = nLevels
    else
       chunklevels = min(l2gp%nLevels, 500)     ! was .., 5)
       chunklevels = max(chunklevels, 1)
    endif
    
    ! Create the swath within the file
    ! print*,"Creating swath called ",name
    swid = HE5_SWcreate(L2FileHandle, trim(name))
    !print*,"Swath ",name,"has SW id :",swid
    if ( swid == -1 ) then
       msr = 'Failed to create swath ' // TRIM(name) &
        & // ' (maybe has the same name as another swath in this file?)'
    end if

    ! Define dimensions

    ! Defining special "unlimited dimension called UNLIM
    ! print*,"Defined Unlim with size", HE5S_UNLIMITED_f
    ! status = HE5_SWdefdim(swid, UNLIM, HE5S_UNLIMITED_F)
    status = mls_swdefdim(swid, UNLIM, HE5S_UNLIMITED_F, &
      & hdfVersion=HDFVERSION_5)

    status = mls_swdefdim(swid, 'nTimes', l2gp%nTimes, &
      & hdfVersion=HDFVERSION_5)

    if ( l2gp%nLevels > 0 ) then
      status = mls_swdefdim(swid, 'nLevels', l2gp%nLevels, &
        & hdfVersion=HDFVERSION_5)
    end if

    if ( l2gp%nFreqs > 0 ) then
      status = mls_swdefdim(swid, 'nFreqs', l2gp%nFreqs, &
        & hdfVersion=HDFVERSION_5)
    end if

    ! Define horizontal geolocation fields using above dimensions

!    print*,"Defining geolocation field ",GEO_FIELD1," of dim. ", DIM_NAME1
!    print*,"... and of type ",HE5T_NATIVE_FLOAT
    chunk_rank=1
    chunk_dims=1
    chunk_dims(1)=CHUNKTIMES
    status = mls_gfldsetup(swid, 'Latitude', 'nTimes', MAX_DIML1, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)
!    print*,"Defined geolocation field ",GEO_FIELD1,"of dim.", DIM_NAME1
!    print*,"... and of type ",HE5T_NATIVE_FLOAT

    status = mls_gfldsetup(swid, 'Longitude', 'nTimes', MAX_DIML1, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)

    status = mls_gfldsetup(swid, 'Time', 'nTimes', MAX_DIML1, &
      & HE5T_NATIVE_DOUBLE, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)

    status = mls_gfldsetup(swid, 'LocalSolarTime', 'nTimes', MAX_DIML1, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)

    status = mls_gfldsetup(swid, 'SolarZenithAngle', 'nTimes', MAX_DIML1, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)

    status = mls_gfldsetup(swid, 'LineOfSightAngle', 'nTimes', MAX_DIML1, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)

    status = mls_gfldsetup(swid, 'OrbitGeodeticAngle', 'nTimes', MAX_DIML1, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)

    status = mls_gfldsetup(swid, 'ChunkNumber', 'nTimes', MAX_DIML1, &
      & HE5T_NATIVE_INT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)

    if ( l2gp%nLevels > 0 ) then

      status = mls_gfldsetup(swid, 'Pressure', 'nLevels', MAX_DIML, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, 0, chunk_dims, &
      & hdfVersion=HDFVERSION_5)
    end if

    if ( l2gp%nFreqs > 0 ) then

      status = mls_gfldsetup(swid, 'Frequency', 'nFreqs', MAX_DIML, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, 0, chunk_dims, &
      & hdfVersion=HDFVERSION_5)
    end if

    ! Define data fields using above dimensions

    if ( (l2gp%nFreqs > 0) .and. (l2gp%nLevels > 0) ) then
       chunk_rank=3
       chunk_dims(1:3)=(/ CHUNKFREQS,CHUNKLEVELS,CHUNKTIMES /)

      status = mls_dfldsetup(swid, 'L2gpValue', 'nFreqs,nLevels,nTimes', &
      & MAX_DIML123, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)

      status = mls_dfldsetup(swid, 'L2gpPrecision', 'nFreqs,nLevels,nTimes', &
      & MAX_DIML123, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)

    else if ( l2gp%nLevels > 0 ) then
       chunk_rank=2
       chunk_dims(1:7)=(/ CHUNKLEVELS,CHUNKTIMES,37,38,39,47,49/)

      status = mls_dfldsetup(swid, 'L2gpValue', 'nLevels,nTimes', &
      & MAX_DIML12, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)

      status = mls_dfldsetup(swid, 'L2gpPrecision', 'nLevels,nTimes', &
      & MAX_DIML12, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)

    else
       chunk_rank=1
       chunk_dims(1)=CHUNKTIMES

      status = mls_dfldsetup(swid, 'L2gpValue', 'nTimes', &
      & MAX_DIML1, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)

      status = mls_dfldsetup(swid, 'L2gpPrecision', 'nTimes', &
      & MAX_DIML1, &
      & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=HDFVERSION_5)

    end if

!    print*,"Defining data field ",DATA_FIELD3,"of dim.", DIM_NAME1
!    print*,"... and of type ",HE5T_NATIVE_CHAR

    chunk_rank=1
    chunk_dims(1)=CHUNKTIMES
    status = mls_dfldsetup(swid, 'Status', 'nTimes', &
    & MAX_DIML1, &
    & HE5T_NATIVE_CHAR, HDFE_NOMERGE, chunk_rank, chunk_dims, &
    & hdfVersion=HDFVERSION_5)
!    status = HE5_SWdefdfld(swid, DATA_FIELD3, DIM_NAME1,MAX_DIML1,&
!         HE5T_NATIVE_CHAR, HDFE_NOMERGE)
!    IF ( status == -1 ) THEN
!       msr = DAT_ERR // DATA_FIELD3
!       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
!    END IF

!    print*,"Defined data field ",DATA_FIELD3,"of dim.", DIM_NAME1
    chunk_rank=1
    chunk_dims(1)=CHUNKTIMES
    status = mls_dfldsetup(swid, 'Quality', 'nTimes', &
    & MAX_DIML1, &
    & HE5T_NATIVE_FLOAT, HDFE_NOMERGE, chunk_rank, chunk_dims, &
    & hdfVersion=HDFVERSION_5)

    ! Set Fill Values
    status =  mls_swsetfill(swid, DATA_FIELDS, HE5T_NATIVE_FLOAT, &
     & l2gp%MissingValue)
    if ( status == -1 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
         & 'Cant set fill for ' // DATA_FIELDS )
    status =  mls_swsetfill(swid, GEO_FIELDS, HE5T_NATIVE_FLOAT, &
     & l2gp%MissingValue)
    if ( status == -1 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
         & 'Cant set fill for ' // GEO_FIELDS )
    status =  mls_swsetfill(swid, 'Time', HE5T_NATIVE_DOUBLE, &
     & real(l2gp%MissingValue, r8) )
    if ( status == -1 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
         & 'Cant set fill for ' // 'Time' )
    status =  mls_swsetfill(swid, 'ChunkNumber', HE5T_NATIVE_INT, &
     & int(l2gp%MissingValue) )
    if ( status == -1 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
         & 'Cant set fill for ' // 'ChunkNumber' )
    if ( l2gp%nLevels > 0 ) then
      status =  mls_swsetfill(swid, 'Pressure', HE5T_NATIVE_FLOAT, &
       & l2gp%MissingValue)
      if ( status == -1 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
           & 'Cant set fill for ' // 'Pressure' )
    endif
    if ( l2gp%nFreqs > 0 ) then
      status =  mls_swsetfill(swid, 'Frequency', HE5T_NATIVE_FLOAT, &
       & l2gp%MissingValue)
      if ( status == -1 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
           & 'Cant set fill for ' // 'Frequency' )
    endif

    ! Detach from the HE5_SWath interface.This stores the swath info within the
    ! file and must be done before writing or reading data to or from the
    ! swath. (May be un-necessary for HDF5 -- test program works OK without.)
    ! 
    status = HE5_SWdetach(swid)
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Failed to detach from swath interface after definition.' )
    end if

    !--------------------------------------
  end subroutine OutputL2GP_createFile_hdf5
  !--------------------------------------

  !-----------------------------------------  OutputL2GP_writeGeo_hdf5  -----
  subroutine OutputL2GP_writeGeo_hdf5 (l2gp, l2FileHandle, swathName,offset)

  use HDFEOS5, only: HE5_swattach, HE5_swdetach
  use MLSHDFEOS, only: mls_swwrfld
    ! Brief description of subroutine
    ! This subroutine writes the geolocation fields to an L2GP output file.

    ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp
    integer, intent(in) :: l2FileHandle ! From swopen
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

    swid = HE5_SWattach (l2FileHandle, name)

    ! Write data to the fields

    stride = 1
    start = myOffset ! Please do not set to zero
    edge(1) = l2gp%nTimes
    status = mls_SWwrfld(swid, 'Latitude', start, stride, edge, &
         real(l2gp%latitude), hdfVersion=HDFVERSION_5)

    status = mls_SWwrfld(swid, 'Longitude', start, stride, edge, &
         real(l2gp%longitude), hdfVersion=HDFVERSION_5)

    status = mls_SWwrfld(swid, 'Time', start, stride, edge, &
         l2gp%time, hdfVersion=HDFVERSION_5)

    status = mls_SWwrfld(swid, 'LocalSolarTime', start, stride, edge, &
        real(l2gp%solarTime), hdfVersion=HDFVERSION_5)

    status = mls_SWwrfld(swid, 'SolarZenithAngle', start, stride, edge, &
         real(l2gp%solarZenith), hdfVersion=HDFVERSION_5)

    status = mls_SWwrfld(swid, 'LineOfSightAngle', start, stride, edge, &
         real(l2gp%losAngle), hdfVersion=HDFVERSION_5)

    status = mls_SWwrfld(swid, 'OrbitGeodeticAngle', start, stride, edge, &
         real(l2gp%geodAngle), hdfVersion=HDFVERSION_5)

    status = mls_SWwrfld(swid, 'ChunkNumber', start, stride, edge, &
         l2gp%chunkNumber, hdfVersion=HDFVERSION_5)

    if ( l2gp%nLevels > 0 ) then
       edge(1) = l2gp%nLevels
       start(1)=0 ! needed because offset may have made this /=0
       status = mls_SWwrfld(swid, 'Pressure', start, stride, edge, &
            real(l2gp%pressures), hdfVersion=HDFVERSION_5)
    end if

    if ( l2gp%nFreqs > 0 ) then
       edge(1) = l2gp%nFreqs
       start(1)=0 ! needed because offset may have made this /=0
       status = mls_SWwrfld(swid, 'Frequency', start, stride, edge, &
            real(l2gp%frequency), hdfVersion=HDFVERSION_5)
    end if

    ! Detach from the swath interface.  

    status = HE5_SWdetach(swid)
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach from swath interface' )
    end if

    !------------------------------------
  end subroutine OutputL2GP_writeGeo_hdf5
  !------------------------------------

  !----------------------------------------  OutputL2GP_writeData_hdf5  -----
  subroutine OutputL2GP_writeData_hdf5(l2gp, l2FileHandle, swathName,offset)

  use HDFEOS5, only: HE5_swattach, HE5_swdetach
  use MLSHDFEOS, only: mls_swwrfld
    ! Brief description of subroutine
    ! This subroutine writes the data fields to an L2GP output file.
    ! For now, you have to write all of l2gp, but you can choose to write
    ! it at some offset into the file
    ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp
    integer, intent(in) :: l2FileHandle ! From swopen
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
    integer, allocatable, dimension(:,:) :: string_buffer

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
    swid = HE5_SWattach (l2FileHandle, name)
    if ( l2gp%nFreqs > 0 ) then
       ! Value and Precision are 3-D fields
       status = mls_SWwrfld(swid, 'L2gpValue', start, stride, edge, &
            & reshape(real(l2gp%l2gpValue), (/size(l2gp%l2gpValue)/)), &
            & hdfVersion=HDFVERSION_5 )
       status = mls_SWwrfld(swid, 'L2gpPrecision', start, stride, edge, &
            & reshape(real(l2gp%l2gpPrecision), (/size(l2gp%l2gpPrecision)/)), &
            & hdfVersion=HDFVERSION_5 )

    else if ( l2gp%nLevels > 0 ) then
       ! Value and Precision are 2-D fields
      
       status = mls_SWwrfld( swid, 'L2gpValue', start(2:3), stride(2:3), &
            edge(2:3), real(l2gp%l2gpValue(1,:,:) ), hdfVersion=HDFVERSION_5)
       status = mls_SWwrfld( swid, 'L2gpPrecision', start(2:3), stride(2:3), &
            edge(2:3), real(l2gp%l2gpPrecision(1,:,:) ), hdfVersion=HDFVERSION_5)
    else

       ! Value and Precision are 1-D fields
       status = mls_SWwrfld( swid, 'L2gpValue', start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpValue(1,1,:) ), hdfVersion=HDFVERSION_5)
       status = mls_SWwrfld( swid, 'L2gpPrecision', start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpPrecision(1,1,:) ), hdfVersion=HDFVERSION_5)
    end if

    ! 1-D status & quality fields

    !HDF-EOS5 won't write a dataset of chars from FORTRAN

   if(USEINTS4STRINGS) then
!      allocate(string_buffer(1,l2gp%nTimes))
!      call strings2Ints(l2gp%status, string_buffer)
!      status = mls_swwrfld(swid,DATA_FIELD3,start(3:3),stride(3:3),edge(3:3),&
!           string_buffer, hdfVersion=HDFVERSION_5)
!      deallocate(string_buffer)
    else
      !    status = HE5_SWwrfld(swid, DATA_FIELD3, start(3:3), stride(3:3),&
      !        edge(3:3), l2gp%status) ! 
      status=0
      call MLSMessage(MLSMSG_Debug, ModuleName, &
        "writing of status field disabled")

      if ( status == -1 ) then
         msr = WR_ERR // DATA_FIELD3
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
    end if
    status = mls_SWwrfld(swid, 'Quality', start(3:3), stride(3:3), edge(3:3), &
         real(l2gp%quality), hdfVersion=HDFVERSION_5)

    !     Detach from the swath interface.

    status = HE5_SWdetach(swid)
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach  from swath interface' )
    end if


    !-------------------------------------
  end subroutine OutputL2GP_writeData_hdf5
  !-------------------------------------

  !----------------------------------------  OutputL2GP_attributes_hdf5  -----
  subroutine OutputL2GP_attributes_hdf5(l2gp, l2FileHandle, swathName)

  use HDFEOS5, only: HE5T_NATIVE_INT, HE5T_NATIVE_REAL, HE5T_NATIVE_DOUBLE, &
    & HE5T_NATIVE_SCHAR, &
    & HE5_SWattach, HE5_SWdetach
  use he5_swapi, only: he5_swwrattr, he5_swwrlattr
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
      & 'Latitude,Longitude,Time,LocalSolarTime,SolarZenithAngle,LineOfSightAngle,OrbitGeodeticAngle,ChunkNumber,Pressure,Frequency'
    character (len=*), parameter :: GeolocationUnits = &
      & 'deg,deg,s,h,deg,deg,deg,NoUnits,hPa,GHz'
    character (len=*), parameter :: GeoUniqueFieldDefinition = &
      & 'HMT,HMT,AS,HMT,HMT,M,M,M,AS,M'   ! These are abbreviated values
    character (len=*), parameter :: UniqueFieldDefKeys = &
      & 'HM,HMT,MT,AS,M'
    character (len=*), parameter :: UniqueFieldDefValues = &
      & 'HIRDLS-MLS-Shared,HIRDLS-MLS-TES-Shared,MLS-TES-Shared,Aura-Shared,MLS-Specific'  ! Expanded values
    character (len=*), parameter :: Species = &
      & 'Temperature,BrO,CH3CN,CO,ClO,GPH,HCl,HCN,H2O,H2O2,HNO3,HOCl,HO2,N2,N2O,OH,O2,O3,RHI,SO2'
    character (len=*), parameter :: SpUniqueFieldDefinition = &
      & 'HMT,M,M,MT,M,M,M,M,HMT,M,HMT,M,M,M,HM,M,M,HMT,M,M'   ! These are abbreviated values

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
    ! swid = HE5_SWattach (l2FileHandle, name)
    !print*," attached swath with swid=",swid," filehandle=",l2FileHandle
    !status = he5_swwrattr(swid, 'Instrument Name', H5T_NATIVE_CHARACTER, &
    !  & 1, 'MLS Aura')
    ! if ( status == -1 ) then
    !   call MLSMessage ( MLSMSG_Warning, ModuleName, &
    !        & 'Failed to write swath attribute' )
    !     Detach from the swath interface.
    ! call sw_writeglobalattr(swid)
    ! print *, 'Writing global attributes'
    call he5_writeglobalattr(l2FileHandle)

    swid = HE5_SWattach (l2FileHandle, name)
    
    !   - -   S w a t h   A t t r i b u t e s   - -
    ! print *, 'Writing swath attributes'
    status = he5_swwrattr(swid, 'Pressure', rgp_type, size(l2gp%pressures), &
      & l2gp%pressures)
    field_name = 'Pressure'
    status = he5_swwrattr(swid, 'VerticalCoordinate', HE5T_NATIVE_SCHAR, 1, &
      & field_name)
    status = he5_swwrattr(swid, 'MissingValue', rgp_type, 1, &
      & (/ real(l2gp%MissingValue, rgp) /) )
    
    !   - -   G e o l o c a t i o n   A t t r i b u t e s   - -
    ! print *, 'Writing geolocation attributes'
    do field=1, NumGeolocFields
      ! Take care not to write attributes to "missing fields"
      if ( trim(theTitles(field)) == 'Frequency' &
        & .and. l2gp%nFreqs < 1 ) then
        field_name = ''
      elseif ( trim(theTitles(field)) == 'Pressure' &
        & .and. l2gp%nLevels < 1 ) then
        field_name = ''
      else
        ! print *, 'field ', field
        ! print *, 'title ', trim(theTitles(field))
        ! print *, 'units ', trim(theUnits(field))
        call GetStringHashElement (GeolocationTitles, &
          & GeoUniqueFieldDefinition, trim(theTitles(field)), &
          & abbr_uniq_fdef, .false.)
        call GetStringHashElement (UniqueFieldDefKeys, &
          & UniqueFieldDefValues, trim(abbr_uniq_fdef), &
          & expnd_uniq_fdef, .false.)
        status = he5_swwrlattr(swid, trim(theTitles(field)), 'Title', &
          & HE5T_NATIVE_SCHAR, 1, theTitles(field))
        status = he5_swwrlattr(swid, trim(theTitles(field)), 'Units', &
          & HE5T_NATIVE_SCHAR, 1, theUnits(field))

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
        status = he5_swwrlattr(swid, trim(theTitles(field)), &
          & 'UniqueFieldDefinition', &
          & HE5T_NATIVE_SCHAR, 1, trim(expnd_uniq_fdef))
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
    ! print *, 'Writing data attributes'
    ! print *, 'title ', trim(field_name)
    ! print *, 'units ', trim(units_name)
    status = he5_swwrlattr(swid, 'L2gpValue', 'Title', &
      & HE5T_NATIVE_SCHAR, 1, field_name)
    status = he5_swwrlattr(swid, 'L2gpValue', 'Units', &
      & HE5T_NATIVE_SCHAR, 1, units_name)
    status = he5_swwrlattr(swid, 'L2gpValue', 'MissingValue', &
      & rgp_type, 1, (/ real(l2gp%MissingValue, rgp) /) )
    status = he5_swwrlattr(swid, 'L2gpValue', &
      & 'UniqueFieldDefinition', &
      & HE5T_NATIVE_SCHAR, 1, trim(expnd_uniq_fdef))
    status = he5_swwrlattr(swid, 'L2gpPrecision', 'Title', &
      & HE5T_NATIVE_SCHAR, 1, trim(field_name)//'Precision')
    status = he5_swwrlattr(swid, 'L2gpPrecision', 'Units', &
      & HE5T_NATIVE_SCHAR, 1, units_name)
    status = he5_swwrlattr(swid, 'L2gpPrecision', 'MissingValue', &
      & rgp_type, 1, (/ real(l2gp%MissingValue, rgp) /) )
    status = he5_swwrlattr(swid, 'L2gpPrecision', &
      & 'UniqueFieldDefinition', &
      & HE5T_NATIVE_SCHAR, 1, trim(expnd_uniq_fdef))

    ! ('Status' data field not yet written)
    
    status = he5_swwrlattr(swid, 'Quality', 'Title', &
      & HE5T_NATIVE_SCHAR, 1, trim(field_name)//'Quality')
    status = he5_swwrlattr(swid, 'Quality', 'Units', &
      & HE5T_NATIVE_SCHAR, 1, units_name)
    status = he5_swwrlattr(swid, 'Quality', 'MissingValue', &
      & rgp_type, 1, (/ real(l2gp%MissingValue, rgp) /) )
    status = he5_swwrlattr(swid, 'Quality', &
      & 'UniqueFieldDefinition', &
      & HE5T_NATIVE_SCHAR, 1, 'MLS-Specific')
    
    status = HE5_SWdetach(swid)
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
    sw_id = he5_swattach(l2FileHandle, trim(name))
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
     & trim(name) // ' Precision')
    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & "Error in setting alias from " // TYPE2PRECISIONNAME // &
          & ' to ' // trim(name) // ' Precision' )
    end if
    returnStatus = he5_SWdetach(sw_id)
    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & "Error in detaching swath for setting alias." )
    end if
  !-------------------------------------
  end subroutine SetL2GP_aliases
  !-------------------------------------


  ! This subroutine is an amalgamation of the last three
  ! Should be renamed CreateAndWriteL2GPData
  subroutine WriteL2GPData(l2gp, l2FileHandle, swathName, hdfVersion)

    ! Arguments

    integer, intent(IN) :: l2FileHandle ! From swopen
    type (L2GPData_T), intent(INOUT) :: l2gp
    character (LEN=*), optional, intent(IN) ::swathName!default->l2gp%swathName
    integer, optional, intent(in) :: hdfVersion
    ! Exectuable code

    ! Local
    integer :: myhdfVersion

    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L2GPDEFAULT_HDFVERSION
    endif

    if (myhdfVersion == HDFVERSION_4) then
      call OutputL2GP_createFile_hdf4 (l2gp, l2FileHandle, swathName)
      call OutputL2GP_writeGeo_hdf4 (l2gp, l2FileHandle, swathName)
      call OutputL2GP_writeData_hdf4 (l2gp, l2FileHandle, swathName)
    elseif (myhdfVersion /= HDFVERSION_5) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Unrecognized hdfVersion passed to WriteL2GPData" )
    else
      ! print *, 'Creating hdfeos5 file'
      call OutputL2GP_createFile_hdf5 (l2gp, l2FileHandle, swathName)
      ! print *, 'Writing geolocation data'
      call OutputL2GP_writeGeo_hdf5 (l2gp, l2FileHandle, swathName)
      ! print *, 'Writing values, precision'
      call OutputL2GP_writeData_hdf5 (l2gp, l2FileHandle, swathName)
      ! print *, 'Writing attributes'
      call OutputL2GP_attributes_hdf5 (l2gp, l2FileHandle, swathName)
      call SetL2GP_aliases (l2gp, l2FileHandle, swathName)
    endif

  end subroutine WriteL2GPData
  !-------------------------------------------------------------


  subroutine AppendL2GPData(l2gp,l2FileHandle,swathName,offset, hdfVersion)
    ! sticks l2gp into the swath swathName in the file pointed at by
    ! l2FileHandle,starting at the profile number "offset" (First profile
    ! in the file has offset==0). If this runs off the end of the swath, 
    ! it is lengthened automagically. 
    ! Arguments

    integer, intent(IN) :: l2FileHandle ! From swopen
    type (L2GPData_T), intent(INOUT) :: l2gp
    character (LEN=*), optional, intent(IN) ::swathName!default->l2gp%swathName
    integer,intent(IN),optional::offset
    integer, optional, intent(in) :: hdfVersion
    ! Local
    integer :: myhdfVersion

    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L2GPDEFAULT_HDFVERSION
    endif

    if (myhdfVersion == HDFVERSION_4) then
      call OutputL2GP_writeGeo_hdf4 (l2gp, l2FileHandle, swathName, offset)
      call OutputL2GP_writeData_hdf4 (l2gp, l2FileHandle, swathName, offset)
    elseif (myhdfVersion /= HDFVERSION_5) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Unrecognized hdfVersion passed to AppendL2GPData" )
    else
      call OutputL2GP_writeGeo_hdf5 (l2gp, l2FileHandle, swathName, offset)
      call OutputL2GP_writeData_hdf5 (l2gp, l2FileHandle, swathName, offset)
    endif

  end subroutine AppendL2GPData

  ! ------------------------------------------ DUMP_L2GP_DATABASE ------------

  subroutine DUMP_L2GP_DATABASE ( L2gp, Name, ColumnsOnly, Details )

    ! Dummy arguments
    type (l2gpData_T), intent(in) ::          L2GP(:)
    character(len=*), intent(in), optional :: Name
    logical, intent(in), optional ::          ColumnsOnly ! if true, dump only with columns
    integer, intent(in), optional :: DETAILS

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
      call dump(l2gp(i), ColumnsOnly, Details)
    end do
      
  end subroutine DUMP_L2GP_DATABASE

  ! ------------------------------------------ DUMP_L2GP ------------


  subroutine Dump_L2GP ( L2gp, ColumnsOnly, Details )

    ! Dummy arguments
    type (l2gpData_T), intent(in) ::          L2GP
    logical, intent(in), optional ::          ColumnsOnly ! if true, dump only with columns
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 1

    ! Local variables
    integer :: ierr
    logical :: myColumnsOnly
    integer :: MYDETAILS
    real(r8) :: FillValue
    integer :: ChunkFillValue

    ! Executable code
    myDetails = 1
    if ( present(details) ) myDetails = details
    
    if( present(ColumnsOnly)) then
      myColumnsOnly = ColumnsOnly
    else
      myColumnsOnly = .false.
    endif

      if ( myColumnsOnly .and. l2gp%nLevels > 1 ) return
      
      FillValue = real(l2gp%MissingValue, r8)
      ChunkFillValue = int(l2gp%MissingValue)

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
      if ( myDetails < -1 ) return
      call output ( 'nTimes: ')
      call output ( l2gp%nTimes, 5)
      call output ( '  nLevels: ')
      call output ( l2gp%nLevels, 3)
      call output ( '  nFreqs: ')
      call output ( l2gp%nFreqs, 3, advance='yes')
      if ( myDetails < 0 ) return
      call dump ( l2gp%pressures, 'Pressures:' )
      
      call dump ( l2gp%latitude, 'Latitude:' )
      
      call dump ( l2gp%longitude, 'Longitude:' )
      
      call dump ( l2gp%solarTime, 'SolarTime:' )
      
      call dump ( l2gp%solarZenith, 'SolarZenith:' )
      
      call dump ( l2gp%losAngle, 'LOSAngle:' )
      
      call dump ( l2gp%geodAngle, 'geodAngle:' )
      
      call dump ( l2gp%time, 'Time:' )
      
      call dump ( l2gp%chunkNumber, 'ChunkNumber:' )
      
      if ( associated(l2gp%frequency) ) &
        & call dump ( l2gp%frequency, 'Frequencies:' )
      
      if ( myDetails < 1 ) return
      call dump ( real(l2gp%l2gpValue, r8), 'L2GPValue:', &
        & FillValue=FillValue )
      
      call dump ( real(l2gp%l2gpPrecision, r8), 'L2GPPrecision:', &
        & FillValue=FillValue )
      
      !    call dump ( l2gp%status, 'Status:' )
      
      call dump ( l2gp%quality, 'Quality:' )
      
  end subroutine Dump_L2GP
    
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
! Revision 2.62  2003/04/15 23:14:50  pwagner
! New chunking so hdfeos5 files dont balloon; some swsetfill tweaking
!
! Revision 2.61  2003/04/11 23:35:10  pwagner
! Added new UniqueFieldDefinition attribute; sets fill and MissingValue attributes for all fields
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
