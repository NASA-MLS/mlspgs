! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module L2GPData                 ! Creation, manipulation and I/O for L2GP Data
!=============================================================================
  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use DUMP_0, only: DUMP
  use Hdf, only: DFNT_CHAR8, DFNT_FLOAT32, DFNT_INT32, DFNT_FLOAT64
  use HDF5_params
  use HDFEOS!, only: SWATTACH, SWCREATE, SWDEFDFLD, SWDEFDIM, SWDEFGFLD, &
     !& SWDETACH
  use HDFEOS5
  use HE5_SWAPI 
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
       & MLSMSG_Error, MLSMSG_Warning
  use MLSCommon, only: R8
  use MLSStrings, only: ints2Strings, strings2Ints
  use OUTPUT_M, only: OUTPUT ! Added as HDF4 version uses it
  use SWAPI, only: SWWRFLD, SWRDFLD
  use STRING_TABLE, only: DISPLAY_STRING

  implicit none

  private
  public :: L2GPData_T
  public :: L2GPNameLen
  public :: AddL2GPToDatabase,  DestroyL2GPContents,  DestroyL2GPDatabase, &
    & Dump, ExpandL2GPDataInPlace, AppendL2GPData, &
    & ReadL2GPData, SetupNewL2GPRecord,  WriteL2GPData !, &
!    & OutputL2GP_createFile, OutputL2GP_writeData,  OutputL2GP_writeGeo
  !INTEGER:: HE5_SWRDFLD
  !EXTERNAL HE5_SWRDFLD !Should USE HE5_SWAPI
  !---------------------------- RCS Ident Info -------------------------------
  character(len=256), private :: Id = &
       & "$Id$"
  character(len=*), parameter, private :: ModuleName = &
       & "$RCSfile$"
  !---------------------------------------------------------------------------

  interface DUMP !And this does WTF? 
    module procedure DUMP_L2GP
    module procedure DUMP_L2GP_DataBase
  end interface

  ! This module defines datatypes and gives basic routines for storing and
  ! manipulating L2GP data.

  ! First some local parameters

  ! Assume L2GP files w/o explicit hdfVersion field are this
  ! 4 corresponds to hdf4, 5 to hdf5 in L2GP, L2AUX, etc. 
  integer, parameter :: L2GPDEFAULT_HDFVERSION = 4            

  integer, parameter :: L2GPNameLen = 80

   character (len=*), parameter :: DATA_FIELD1 = 'L2gpValue'
   character (len=*), parameter :: DATA_FIELD2 = 'L2gpPrecision'
   character (len=*), parameter :: DATA_FIELD3 = 'Status'
   character (len=*), parameter :: DATA_FIELD4 = 'Quality'

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

   character (len=*), parameter :: DIM_NAME1 = 'nTimes'
   character (len=*), parameter :: DIM_NAME2 = 'nLevels'
   character (len=*), parameter :: DIM_NAME3 = 'nFreqs'
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

     ! Now the dimensions of the data

     integer :: nTimes          ! Total number of profiles
     integer :: nLevels         ! Total number of surfaces (==1 for col. abund)
     integer :: nFreqs          ! Number of frequencies in breakdown

     ! Now we store the geolocation fields, first the vertical one:
     ! (The following has the tropopause if the swath is a column abundance)
     real (r8), pointer, dimension(:) :: pressures=>NULL() ! Vertical coords (nLevels)

     ! Now the horizontal geolocation information. Dimensioned (nTimes)
     real (r8), pointer, dimension(:) :: latitude => NULL()
     real (r8), pointer, dimension(:) :: longitude => NULL()
     real (r8), pointer, dimension(:) :: solarTime => NULL()
     real (r8), pointer, dimension(:) :: solarZenith => NULL()
     real (r8), pointer, dimension(:) :: losAngle => NULL()
     real (r8), pointer, dimension(:) :: geodAngle => NULL()
     real (r8), pointer, dimension(:) :: time => NULL()

     integer, pointer, dimension(:) :: chunkNumber=>NULL()

     ! Now we store the `frequency' geolocation field

     real (r8), pointer, dimension(:) :: frequency=>NULL()
     !        dimensioned (nFreqs)

     ! Finally we store the data fields

     real (r8), pointer, dimension(:,:,:) :: l2gpValue=>NULL()
     real (r8), pointer, dimension(:,:,:) :: l2gpPrecision=>NULL()
     ! dimensioned (nFreqs, nLevels, nTimes)

     character (len=1), pointer, dimension(:) :: status=>NULL()
     !                (status is a reserved word in F90)
     real (r8), pointer, dimension(:) :: quality=>NULL()
     ! Both the above dimensioned (nTimes)

  end type L2GPData_T

contains ! =====     Public Procedures     =============================

  !------------------------------------------  SetupNewL2GPRecord  -----
  subroutine SetupNewL2GPRecord ( l2gp, nFreqs, nLevels, nTimes)

    ! This routine sets up the arrays for an l2gp datatype.

    ! Dummy arguments
    type (L2GPData_T), intent(out)  :: l2gp
    integer, intent(in), optional :: nFreqs            ! Dimensions
    integer, intent(in), optional :: nLevels           ! Dimensions
    integer, intent(in), optional :: nTimes            ! Dimensions

    ! Local variables
    integer :: freqsArrayLen, status, surfsArrayLen
    integer :: useNFreqs, useNLevels, useNTimes

    if (present(nFreqs)) then
       useNFreqs=nFreqs
    else
       useNFreqs=0
    endif

    if (present(nLevels)) then
       useNLevels=nLevels
    else
       useNLevels=0
    endif

    if (present(nTimes)) then
       useNTimes=nTimes
    else
       useNTimes=0              ! Default to empty l2gp
    endif

    ! Store the dimensionality

    l2gp%nTimes = useNTimes
    l2gp%nLevels = useNLevels
    l2gp%nFreqs = useNFreqs

    ! But allocate to at least one for times, freqs

    useNTimes=max(useNTimes,1)
    useNLevels=max(useNLevels,1)
    useNFreqs=max(useNFreqs,1)    

    ! Allocate the frequency coordinate

    call allocate_test ( l2gp%pressures, useNLevels, "l2gp%pressures", &
         & ModuleName )

    ! Allocate the vertical coordinate

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

  end subroutine SetupNewL2GPRecord

  !-----------------------------------------  DestroyL2GPContents  -----
  subroutine DestroyL2GPContents ( L2GP )

    ! This routine deallocates all the arrays allocated above.

    ! Dummy arguments
    type (L2GPData_T), intent(inout) :: L2GP
    ! Local variables

    integer status

    ! Executable code

    call deallocate_test ( l2gp%pressures,    "l2gp%pressures",    ModuleName )
    call deallocate_test ( l2gp%latitude,     "l2gp%latitude",     ModuleName )
    call deallocate_test ( l2gp%longitude,    "l2gp%longitude",    ModuleName )
    call deallocate_test ( l2gp%solarTime,    "l2gp%solarTime",    ModuleName )
    call deallocate_test ( l2gp%solarZenith,  "l2gp%solarZenith",  ModuleName )
    call deallocate_test ( l2gp%losAngle,     "l2gp%losAngle",     ModuleName )
    call deallocate_test ( l2gp%losAngle,     "l2gp%losAngle",     ModuleName )
    call deallocate_test ( l2gp%geodAngle,    "l2gp%geodAngle",    ModuleName )
    call deallocate_test ( l2gp%chunkNumber,  "l2gp%chunkNumber",  ModuleName )
    call deallocate_test ( l2gp%time,         "l2gp%time",         ModuleName )
    call deallocate_test ( l2gp%frequency,    "l2gp%frequency",    ModuleName )
    call deallocate_test ( l2gp%l2gpValue,    "l2gp%l2gpValue",    ModuleName )
    call deallocate_test ( l2gp%l2gpPrecision,"l2gp%l2gpPrecision",ModuleName )
    call deallocate_test ( l2gp%status,       "l2gp%status",       ModuleName )
    call deallocate_test ( l2gp%quality,      "l2gp%quality",      ModuleName )
    l2gp%nTimes = 0
    l2gp%nLevels = 0
    l2gp%nFreqs = 0
  end subroutine DestroyL2GPContents

  !---------------------------------------  ExpandL2GPDataInPlace  -----
  subroutine ExpandL2GPDataInPlace ( l2gp, newNTimes )

    ! This subroutine expands an L2GPData_T in place allowing the user to
    ! (1) add more profiles to it; or [or? or WTF?]

    ! Dummy arguments
    type (L2GPData_T), intent(inout) :: l2gp
    integer, optional, intent(in) :: newNTimes

    ! Local variables
    integer :: status                   ! From ALLOCATE
    type (L2GPData_T) :: tempL2gp       ! For copying data around
    integer :: myNTimes, myNColms
    ! Executable code

   if(present(newNTimes)) then
      myNTimes = newNTimes
   else
      myNTimes = l2gp%nTimes
   endif

    ! First do a sanity check


    if ( myNTimes<l2gp%nTimes ) &
         & call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "The number of profiles requested is fewer than those already present" )

    tempL2gp = l2gp ! Copy the pointers to the old information

    ! Now recreate l2gp with the new size.
    ! First, nullify all of the pointers in l2gp, so that a deallocate_test
    ! won't delete them.  After all, we just went to the trouble to preserve
    ! them in TempL2GP!
    nullify ( l2gp%pressures, l2gp%latitude, l2gp%longitude, l2gp%solarTime, &
      & l2gp%solarZenith, l2gp%losAngle, l2gp%losAngle, l2gp%geodAngle, &
      & l2gp%chunkNumber, l2gp%time, l2gp%frequency, l2gp%l2gpValue, &
      & l2gp%l2gpPrecision, l2gp%status, l2gp%quality )
    call SetupNewL2GPRecord( l2gp, nFreqs=l2gp%nFreqs, nLevels=l2gp%nLevels, &
      & nTimes=myNTimes)

    ! Don't forget the `global' stuff
    l2gp%pressures=templ2gp%pressures

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
    type (L2GPData_T), dimension(:), pointer :: tempDatabase
    !This include causes real trouble if you are compiling in a different 
    !directory.
    include "addItemToDatabase.f9h" 

    AddL2GPToDatabase = newSize
  end function AddL2GPToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template database

  subroutine DestroyL2GPDatabase ( DATABASE )

    ! Dummy argument
    type (L2GPData_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer :: l2gpIndex, status

    if ( associated(database)) then
       do l2gpIndex = 1, size(database)
          call DestroyL2GPContents ( database(l2gpIndex) )
       end do
       deallocate ( database, stat=status )
       if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & MLSMSG_deallocate // "database" )
    end if
  end subroutine DestroyL2GPDatabase

  ! ---------------------- ReadL2GPData  -----------------------------

  subroutine ReadL2GPData(L2FileHandle, swathname, l2gp, numProfs, &
       firstProf, lastProf, hdfVersion)
    !------------------------------------------------------------------------

    ! This routine reads an L2GP file, in either hdfVersion,
    ! returning a filled data structure and the !
    ! number of profiles read.

    ! Arguments

    character (len=*), intent(in) :: swathname ! Name of swath
    integer, intent(in) :: L2FileHandle ! Returned by swopen
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( l2GPData_T ), intent(out) :: l2gp ! Result
    integer, intent(out), optional :: numProfs ! Number actually read
    integer, optional, intent(in) :: hdfVersion

    ! Local
    integer :: myhdfVersion

    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L2GPDEFAULT_HDFVERSION
    endif

    if (myhdfVersion == 4) then
      call ReadL2GPData_hdf4(L2FileHandle, swathname, l2gp, numProfs, &
       firstProf, lastProf)
    else
      call ReadL2GPData_hdf5(L2FileHandle, swathname, l2gp, numProfs, &
       firstProf, lastProf)
!       call MLSMessage ( MLSMSG_Error, ModuleName, &
!            & 'This version of L2GPData not yet ready for hdf5' )
    endif

  end subroutine ReadL2GPData

  ! ---------------------- ReadL2GPData_hdf4  -----------------------------

  subroutine ReadL2GPData_hdf4(L2FileHandle, swathname, l2gp, numProfs, &
       firstProf, lastProf)
    !------------------------------------------------------------------------

    ! This routine reads an L2GP file, assuming it is hdfeos2 format
    ! (i.e. hdf4) returning a filled data structure and the
    ! number of profiles read.

    ! Arguments

    character (len=*), intent(in) :: swathname ! Name of swath
    integer, intent(in) :: L2FileHandle ! Returned by swopen
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( l2GPData_T ), intent(out) :: l2gp ! Result
    integer, intent(out), optional :: numProfs ! Number actually read

    ! Local Parameters
    character (len=*), parameter :: SZ_ERR = 'Failed to get size of &
         &dimension '
    character (len=*), parameter :: MLSMSG_INPUT = 'Error in input argument '
    character (len=*), parameter :: MLSMSG_L2GPRead = 'Unable to read L2GP &
                                                     &field:'

    ! Local Variables
    character (len=80) :: list
    character (len=480) :: msr

    integer :: alloc_err, first, freq, lev, nDims, size, swid, status
    integer :: start(3), stride(3), edge(3), dims(3)
    integer :: nFreqs, nLevels, nTimes, nFreqsOr1, nLevelsOr1, myNumProfs
    integer :: nColumns
    integer :: col_start(2), col_stride(2), col_edge(2)

    logical :: firstCheck, lastCheck

    real, allocatable :: realFreq(:), realSurf(:), realProf(:), real3(:,:,:)
!  How to deal with status and columnTypes? swrfld fails
!  with char data on Linux
!  Have recourse to ints2Strings and strings2Ints if USEINTS4STRINGS
!    character (LEN=8), allocatable :: the_status_buffer(:)
!    character (LEN=L2GPNameLen), allocatable :: the_status_buffer(:)
    integer, allocatable, dimension(:,:) :: string_buffer

    ! Attach to the swath for reading

    l2gp%Name = swathname

    swid = swattach(L2FileHandle, TRIM(l2gp%Name))
    if ( swid == -1) call MLSMessage ( MLSMSG_Error, ModuleName, 'Failed to &
         &attach to swath interface for reading.')

    ! Get dimension information

    lev = 0
    freq = 0

    nDims = swinqdims(swid, list, dims)
    if ( nDims == -1) call MLSMessage ( MLSMSG_Error, ModuleName, 'Failed to &
         &get dimension information.')
    if ( INDEX(list,'nLevels') /= 0 ) lev = 1
    if ( INDEX(list,'Freq') /= 0 ) freq = 1

    size = swdiminfo(swid, DIM_NAME1)
    if ( size == -1 ) then
       msr = SZ_ERR // DIM_NAME1
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    l2gp%nTimes = size
    nTimes=size
    if ( lev == 0 ) then
       nLevels = 0
    else
       size = swdiminfo(swid, DIM_NAME2)
       if ( size == -1 ) then
          msr = SZ_ERR // DIM_NAME2
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       nLevels = size

    end if

    if ( freq == 1 ) then
       size = swdiminfo(swid, DIM_NAME3)
       if ( size == -1 ) then
          msr = SZ_ERR // DIM_NAME3
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       nFreqs = size
    else
       nFreqs = 0
    end if

    ! Check optional input arguments

    firstCheck = present(firstProf)
    lastCheck = present(lastProf)

    if ( firstCheck ) then

       if ( (firstProf >= l2gp%nTimes) .OR. (firstProf < 0) ) then
          msr = MLSMSG_INPUT // 'firstProf'
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       else
          first = firstProf
       end if

    else

       first = 0

    end if

    if ( lastCheck ) then

       if ( lastProf < first ) then
          msr = MLSMSG_INPUT // 'lastProf'
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
    & nTimes=myNumProfs )

    ! Allocate temporary arrays

    nFreqsOr1=MAX(nFreqs,1)
    nLevelsOr1=MAX(nLevels, 1)
    allocate ( realProf(myNumProfs), realSurf(l2gp%nLevels), &
      &   string_buffer(1,myNumProfs), &
      &   realFreq(l2gp%nFreqs), &
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

    status = swrdfld(swid, GEO_FIELD1, start(3:3), stride(3:3), edge(3:3), &
      &    realProf)
    if ( status == -1 ) then
       msr = MLSMSG_L2GPRead // GEO_FIELD1
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    l2gp%latitude = realProf

    status = swrdfld(swid, GEO_FIELD2, start(3:3), stride(3:3), edge(3:3), &
      &    realProf)
    if ( status == -1 ) then
       msr = MLSMSG_L2GPRead // GEO_FIELD2
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    l2gp%longitude = realProf

    status = swrdfld(swid, GEO_FIELD3, start(3:3), stride(3:3), edge(3:3), &
      &    l2gp%time)
    if ( status == -1 ) then
       msr = MLSMSG_L2GPRead // GEO_FIELD3
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swrdfld(swid, GEO_FIELD4, start(3:3), stride(3:3), edge(3:3), &
      &    realProf)
    if ( status == -1 ) then
       msr = MLSMSG_L2GPRead // GEO_FIELD4
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    l2gp%solarTime = realProf

    status = swrdfld(swid, GEO_FIELD5, start(3:3), stride(3:3), edge(3:3), &
      &    realProf)
    if ( status == -1 ) then
       msr = MLSMSG_L2GPRead // GEO_FIELD5
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    l2gp%solarZenith = realProf

    status = swrdfld(swid, GEO_FIELD6, start(3:3), stride(3:3), edge(3:3), &
      &    realProf)
    if ( status == -1 ) then
       msr = MLSMSG_L2GPRead // GEO_FIELD6
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    l2gp%losAngle = realProf

    status = swrdfld(swid, GEO_FIELD7, start(3:3), stride(3:3), edge(3:3), &
      &    realProf)
    if ( status == -1 ) then
       msr = MLSMSG_L2GPRead // GEO_FIELD7
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    l2gp%geodAngle = realProf

    status = swrdfld(swid, GEO_FIELD8, start(3:3), stride(3:3), edge(3:3), &
      &    l2gp%chunkNumber)
    if ( status == -1 ) then
       msr = MLSMSG_L2GPRead // GEO_FIELD8
       call MLSMessage ( MLSMSG_Warning, ModuleName, msr )
    end if

    ! Read the pressures vertical geolocation field, if it exists

    if ( lev /= 0 ) then

       status = swrdfld(swid, GEO_FIELD9, start(2:2), stride(2:2), edge(2:2), &
         & realSurf)
       if ( status == -1 ) then
          msr = MLSMSG_L2GPRead // GEO_FIELD9
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

       l2gp%pressures = realSurf

    end if

    ! Read the frequency geolocation field, if it exists

    if ( freq == 1 ) then

       edge(1) = l2gp%nFreqs

       status = swrdfld(swid, GEO_FIELD10, start(1:1), stride(1:1), edge(1:1), &
         & realFreq)
       if ( status == -1 ) then
          msr = MLSMSG_L2GPRead // GEO_FIELD10
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       l2gp%frequency = realFreq

    end if

    ! Read the data fields that may have 1-3 dimensions

    if ( freq == 1 ) then

       status = swrdfld(swid, DATA_FIELD1, start, stride, edge, real3)
       if ( status == -1 ) then
          msr = MLSMSG_L2GPRead // DATA_FIELD1
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       l2gp%l2gpValue = real3

       status = swrdfld(swid, DATA_FIELD2, start, stride, edge, real3)
       if ( status == -1 ) then
          msr = MLSMSG_L2GPRead // DATA_FIELD2
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       l2gp%l2gpPrecision = real3

    else if ( lev == 1 ) then

      status = swrdfld( swid, DATA_FIELD1, start(2:3), stride(2:3), &
        &   edge(2:3), real3(1,:,:) )
      if ( status == -1 ) then
        msr = MLSMSG_L2GPRead // DATA_FIELD1
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
      l2gp%l2gpValue = real3
      
      status = swrdfld( swid, DATA_FIELD2, start(2:3), stride(2:3), &
        & edge(2:3), real3(1,:,:) )
      if ( status == -1 ) then
        msr = MLSMSG_L2GPRead // DATA_FIELD2
        call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
      l2gp%l2gpPrecision = real3
      
    else

       status = swrdfld( swid, DATA_FIELD1, start(3:3), stride(3:3), edge(3:3), &
         &   real3(1,1,:) )
       if ( status == -1 ) then
          msr = MLSMSG_L2GPRead // DATA_FIELD1
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       l2gp%l2gpValue = real3

       status = swrdfld( swid, DATA_FIELD2, start(3:3), stride(3:3), edge(3:3), &
            real3(1,1,:) )
       if ( status == -1 ) then
          msr = MLSMSG_L2GPRead // DATA_FIELD2
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
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


    l2gp%status = ' ' ! So it has a value.

   if(USEINTS4STRINGS) then
       status = swrdfld(swid, DATA_FIELD3,start(3:3),stride(3:3),edge(3:3),&
         string_buffer)
       if ( status == -1 ) then
         msr = MLSMSG_L2GPRead // DATA_FIELD3
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       call ints2Strings(string_buffer, l2gp%status)
    end if

    status = swrdfld(swid, DATA_FIELD4, start(3:3), stride(3:3), edge(3:3), &
      &   realProf)
    if ( status == -1 ) then
       msr = MLSMSG_L2GPRead // DATA_FIELD4
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
    l2gp%quality = realProf

    ! Deallocate local variables

    deallocate ( realFreq, realSurf, realProf, real3, STAT=alloc_err )
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

  subroutine ReadL2GPData_hdf5(L2FileHandle, swathname, l2gp, numProfs, &
       firstProf, lastProf)
    !------------------------------------------------------------------------

    ! This routine reads an L2GP file, returning a filled data structure and the !
    ! number of profiles read.

    ! Arguments

    character (LEN=*), intent(IN) :: swathname ! Name of swath
    integer, intent(IN) :: L2FileHandle ! Returned by swopen
    integer, intent(IN), optional :: firstProf, lastProf ! Defaults to first and last
    type( L2GPData_T ), intent(OUT) :: l2gp ! Result
    integer, intent(OUT),optional :: numProfs ! Number actually read

    ! Local Parameters
    character (LEN=*), parameter :: SZ_ERR = 'Failed to get size of &
         &dimension '
    character (LEN=*), parameter :: MLSMSG_INPUT = 'Error in input argument '
    character (LEN=*), parameter :: MLSMSG_L2GPRead = 'Unable to read L2GP &
                                                     &field:'
    ! Local Variables
    character (LEN=80) :: list
    character (LEN=480) :: msr

    integer :: alloc_err, first, freq, lev, nDims, size, swid, status
    integer :: start(3), stride(3), edge(3), dims(3)
    integer :: nFreqs, nLevels, nTimes, nFreqsOr1, nLevelsOr1, myNumProfs
    integer :: nColumns
    integer :: col_start(2), col_stride(2), col_edge(2)

    logical :: firstCheck, lastCheck

    real, allocatable :: realFreq(:), realSurf(:), realProf(:), real3(:,:,:)
!  How to deal with status and columnTypes? swrfld fails
!  with char data on Linux with HDF4. With HDF5 we may or may not need to
!  Have recourse to ints2Strings and strings2Ints if USEINTS4STRINGS
!    character (LEN=8), allocatable :: the_status_buffer(:)
!    character (LEN=L2GPNameLen), allocatable :: the_status_buffer(:)
    integer, allocatable, dimension(:,:) :: string_buffer

    ! Attach to the swath for reading

    l2gp%Name = swathname
    
    swid = HE5_SWattach(L2FileHandle, l2gp%Name)
    if (swid == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &attach to swath interface for reading.')

    ! Get dimension information

    lev = 0
    freq = 0

    nDims = HE5_SWinqdims(swid, list, dims)
    if (nDims == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &get dimension information.')
    if ( index(list,'nLevels') /= 0 ) lev = 1
    if ( index(list,'Freq') /= 0 ) freq = 1

    size = HE5_SWdiminfo(swid, DIM_NAME1)
    if (size == -1) then
       msr = SZ_ERR // DIM_NAME1
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif
    l2gp%nTimes = size
    nTimes=size
    if (lev == 0) then
       nLevels = 0
    else
       size = HE5_SWdiminfo(swid, DIM_NAME2)
       if (size == -1) then
          msr = SZ_ERR // DIM_NAME2
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif
       nLevels = size

    endif

    if (freq == 1) then
       size = HE5_SWdiminfo(swid, DIM_NAME3)
       if (size == -1) then
          msr = SZ_ERR // DIM_NAME3
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif
       nFreqs = size
    else
       nFreqs = 0
    endif

    ! Check optional input arguments

    firstCheck = present(firstProf)
    lastCheck = present(lastProf)

    if (firstCheck) then

       if ( (firstProf >= l2gp%nTimes) .or. (firstProf < 0) ) then
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
      &  nTimes=mynumProfs)

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

    status = HE5_SWrdfld(swid, GEO_FIELD1, start(3:3), stride(3:3), &
      edge(3:3), realProf)
    if (status == -1) then
       msr = MLSMSG_L2GPRead // GEO_FIELD1
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif
    l2gp%latitude = realProf

    status = HE5_SWrdfld(swid, GEO_FIELD2, start(3:3), stride(3:3), edge(3:3),&
      &    realProf)
    if (status == -1) then
       msr = MLSMSG_L2GPRead // GEO_FIELD2
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif
    l2gp%longitude = realProf

    status = HE5_SWrdfld(swid, GEO_FIELD3, start(3:3), stride(3:3), edge(3:3),&
      &    l2gp%time)
    if (status == -1) then
       msr = MLSMSG_L2GPRead // GEO_FIELD3
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif

    status = HE5_SWrdfld(swid, GEO_FIELD4, start(3:3), stride(3:3), edge(3:3),&
      &    realProf)
    if (status == -1) then
       msr = MLSMSG_L2GPRead // GEO_FIELD4
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif
    l2gp%solarTime = realProf

    status = HE5_SWrdfld(swid, GEO_FIELD5, start(3:3), stride(3:3), edge(3:3),&
      &    realProf)
    if (status == -1) then
       msr = MLSMSG_L2GPRead // GEO_FIELD5
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif
    l2gp%solarZenith = realProf

    status = HE5_SWrdfld(swid, GEO_FIELD6, start(3:3), stride(3:3), edge(3:3),&
      &    realProf)
    if (status == -1) then
       msr = MLSMSG_L2GPRead // GEO_FIELD6
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif
    l2gp%losAngle = realProf

    status = HE5_SWrdfld(swid, GEO_FIELD7, start(3:3), stride(3:3), edge(3:3),&
      &   realProf)
    if (status == -1) then
       msr = MLSMSG_L2GPRead // GEO_FIELD7
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif
    l2gp%geodAngle = realProf

    status = HE5_SWrdfld(swid, GEO_FIELD8, start(3:3), stride(3:3), edge(3:3),&
      &    l2gp%chunkNumber)
    if (status == -1) then
       msr = MLSMSG_L2GPRead // GEO_FIELD8
       call MLSMessage(MLSMSG_Warning, ModuleName, msr)
    endif

    ! Read the pressures vertical geolocation field, if it exists

    if (lev /= 0) then

       status = HE5_SWrdfld(swid,GEO_FIELD9,start(2:2),stride(2:2), edge(2:2),&
         & realSurf)
       if (status == -1) then
          msr = MLSMSG_L2GPRead // GEO_FIELD9
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif

       l2gp%pressures = realSurf

    endif

    ! Read the frequency geolocation field, if it exists

    if (freq == 1) then

       edge(1) = l2gp%nFreqs

       status = HE5_SWrdfld(swid,GEO_FIELD10,start(1:1),stride(1:1),edge(1:1),&
         & realFreq)
       if (status == -1) then
          msr = MLSMSG_L2GPRead // GEO_FIELD10
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif
       l2gp%frequency = realFreq

    endif

    ! Read the data fields that may have 1-3 dimensions

    if ( freq == 1) then

       status = HE5_SWrdfld(swid, DATA_FIELD1, start, stride, edge, real3)
       if (status == -1) then
          msr = MLSMSG_L2GPRead // DATA_FIELD1
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif
       l2gp%l2gpValue = real3

       status = HE5_SWrdfld(swid, DATA_FIELD2, start, stride, edge, real3)
       if (status == -1) then
          msr = MLSMSG_L2GPRead // DATA_FIELD2
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif
       l2gp%l2gpPrecision = real3

    else if ( lev == 1) then

       status = HE5_SWrdfld( swid, DATA_FIELD1, start(2:3), stride(2:3), &
            edge(2:3), real3 )
!            edge(2:3), real3(1,:,:) )
       if (status == -1) then
          msr = MLSMSG_L2GPRead // DATA_FIELD1
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif
       l2gp%l2gpValue = real3

       status = HE5_SWrdfld( swid, DATA_FIELD2, start(2:3), stride(2:3), &
            edge(2:3), real3(1,:,:) )
       if (status == -1) then
          msr = MLSMSG_L2GPRead // DATA_FIELD2
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif
       l2gp%l2gpPrecision = real3

    else

       status = HE5_SWrdfld(swid,DATA_FIELD1,start(3:3),stride(3:3),edge(3:3),&
         &   real3(1,1,:) )
       if (status == -1) then
          msr = MLSMSG_L2GPRead // DATA_FIELD1
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif
       l2gp%l2gpValue = real3

       status = HE5_SWrdfld( swid, DATA_FIELD2, start(3:3), stride(3:3), edge(3:3), &
            real3(1,1,:) )
       if (status == -1) then
          msr = MLSMSG_L2GPRead // DATA_FIELD2
          call MLSMessage(MLSMSG_Error, ModuleName, msr)
       endif
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

    l2gp%status = ' ' ! So it has a value.

   if(USEINTS4STRINGS) then
       status = HE5_swrdfld(swid,DATA_FIELD3,start(3:3),stride(3:3),edge(3:3),&
         string_buffer)
       if ( status == -1 ) then
         msr = MLSMSG_L2GPRead // DATA_FIELD3
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       call ints2Strings(string_buffer, l2gp%status)
    else
      print*,"Warning: reading of status field disabled"
      status=0
    end if

    status = HE5_SWrdfld(swid, DATA_FIELD4, start(3:3), stride(3:3),&
      edge(3:3),realProf)
    if (status == -1) then
       msr = MLSMSG_L2GPRead // DATA_FIELD4
       call MLSMessage(MLSMSG_Error, ModuleName, msr)
    endif
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

    ! Set numProfs if wanted
    if (present(numProfs)) numProfs=myNumProfs

    !-----------------------------
  end subroutine ReadL2GPData_hdf5
  !-----------------------------

  ! --------------------------------------  OutputL2GP_createFile_hdf4  -----
  subroutine OutputL2GP_createFile_hdf4 (l2gp, L2FileHandle, swathName)

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
       msr = 'Failed to create swath ' // TRIM(name)
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    ! Define dimensions

    status = swdefdim(swid, DIM_NAME1, l2gp%nTimes)
    if ( status == -1 ) then
       msr = DIM_ERR // DIM_NAME1
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    if ( l2gp%nLevels > 0 ) then
       status = swdefdim(swid, DIM_NAME2, l2gp%nLevels)
       if ( status == -1 ) then
          msr = DIM_ERR // DIM_NAME2
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    end if

    if ( l2gp%nFreqs > 1 ) then
       status = swdefdim(swid, DIM_NAME3, l2gp%nFreqs)
       if ( status == -1 ) then
          msr = DIM_ERR // DIM_NAME3
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    end if

    ! Define horizontal geolocation fields using above dimensions

    status = swdefgfld(swid, GEO_FIELD1, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD1
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD2, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD2
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD3, DIM_NAME1, DFNT_FLOAT64, &
         HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD3
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD4, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD4
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD5, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD5
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD6, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD6
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD7, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD7
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefgfld(swid, GEO_FIELD8, DIM_NAME1, DFNT_INT32, &
         HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD8
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    if ( l2gp%nLevels > 0 ) then
       status = swdefgfld(swid, GEO_FIELD9, DIM_NAME2, DFNT_FLOAT32, &
            HDFE_NOMERGE)
       if ( status == -1 ) then
          msr = GEO_ERR // GEO_FIELD9
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    end if

    if ( l2gp%nFreqs > 0 ) then
       status = swdefgfld(swid, GEO_FIELD10, DIM_NAME3, DFNT_FLOAT32, &
            HDFE_NOMERGE)
       if ( status == -1 ) then
          msr = GEO_ERR // GEO_FIELD10
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    end if

    ! Define data fields using above dimensions

    if ( (l2gp%nFreqs > 0) .AND. (l2gp%nLevels > 0) ) then

       status = swdefdfld(swid, DATA_FIELD1, DIM_NAME123, DFNT_FLOAT32, &
            HDFE_NOMERGE)

       if ( status == -1 ) then
          msr = DAT_ERR // DATA_FIELD1 // ' for 3D quantity.'
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if


       status = swdefdfld(swid, DATA_FIELD2, DIM_NAME123, DFNT_FLOAT32, &
            HDFE_NOMERGE)

       if ( status == -1 ) then
          msr = DAT_ERR // DATA_FIELD2 // ' for 3D quantity.'
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if


    else if ( l2gp%nLevels > 0 ) then

       status = swdefdfld(swid, DATA_FIELD1, DIM_NAME12, DFNT_FLOAT32, &
            HDFE_NOMERGE)

       if ( status == -1 ) then
          msr = DAT_ERR // DATA_FIELD1 //  ' for 2D quantity.'
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

       status = swdefdfld(swid, DATA_FIELD2, DIM_NAME12, DFNT_FLOAT32, &
            HDFE_NOMERGE)

       if ( status == -1 ) then
          msr = DAT_ERR // DATA_FIELD2 //  ' for 2D quantity.'
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

    else

       status = swdefdfld(swid, DATA_FIELD1, DIM_NAME1, DFNT_FLOAT32, &
            HDFE_NOMERGE)

       if ( status == -1 ) then
          msr = DAT_ERR // DATA_FIELD1 // ' for 1D quantity.'
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

       status = swdefdfld(swid, DATA_FIELD2, DIM_NAME1, DFNT_FLOAT32, &
            HDFE_NOMERGE)

       if ( status == -1 ) then
          msr = DAT_ERR // DATA_FIELD2 // ' for 1D quantity.'
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

    end if

    status = swdefdfld(swid, DATA_FIELD3, DIM_NAME1, DFNT_CHAR8, &
         HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = DAT_ERR // DATA_FIELD3
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swdefdfld(swid, DATA_FIELD4, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = DAT_ERR // DATA_FIELD4
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

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
  subroutine OutputL2GP_writeGeo_hdf4 (l2gp, l2FileHandle, swathName)

    ! Brief description of subroutine
    ! This subroutine writes the geolocation fields to an L2GP output file.

    ! Arguments

    type( l2GPData_T ), intent(inout) :: l2gp
    integer, intent(in) :: l2FileHandle ! From swopen
    character (len=*), intent(in), optional :: swathName ! Defaults to l2gp%name

    ! Parameters

    character (len=*), parameter :: WR_ERR = &
         & 'Failed to write geolocation field '

    ! Variables

    character (len=480) :: msr
    character (len=132) :: name ! Either swathName or l2gp%name

    integer :: status, swid
    integer :: start(2), stride(2), edge(2)

    if ( present(swathName) ) then
       name=swathName
    else
       name=l2gp%name
    end if

    swid = swattach (l2FileHandle, name)

    ! Write data to the fields

    stride(1) = 1
    start(1) = 0
    edge(1) = l2gp%nTimes

    status = swwrfld(swid, GEO_FIELD1, start, stride, edge, &
         real(l2gp%latitude))
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD1
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swwrfld(swid, GEO_FIELD2, start, stride, edge, &
         real(l2gp%longitude))
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD2
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swwrfld(swid, GEO_FIELD3, start, stride, edge, &
         l2gp%time)
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD3
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swwrfld(swid, GEO_FIELD4, start, stride, edge, &
        real(l2gp%solarTime))
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD4
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swwrfld(swid, GEO_FIELD5, start, stride, edge, &
         real(l2gp%solarZenith))
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD5
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swwrfld(swid, GEO_FIELD6, start, stride, edge, &
         real(l2gp%losAngle))
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD6
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swwrfld(swid, GEO_FIELD7, start, stride, edge, &
         real(l2gp%geodAngle))
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD7
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = swwrfld(swid, GEO_FIELD8, start, stride, edge, &
         l2gp%chunkNumber)
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD8
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    if ( l2gp%nLevels > 0 ) then
       edge(1) = l2gp%nLevels
       status = swwrfld(swid, GEO_FIELD9, start, stride, edge, &
            real(l2gp%pressures))
       if ( status == -1 ) then
          msr = WR_ERR // GEO_FIELD9
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    end if

    if ( l2gp%nFreqs > 0 ) then
       edge(1) = l2gp%nFreqs
       l2gp%frequency = 0
       status = swwrfld(swid, GEO_FIELD10, start, stride, edge, &
            real(l2gp%frequency))
       if ( status == -1 ) then
          msr = WR_ERR // GEO_FIELD10
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
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
  subroutine OutputL2GP_writeData_hdf4(l2gp, l2FileHandle, swathName)

    ! Brief description of subroutine
    ! This subroutine writes the data fields to an L2GP output file.

    ! Arguments

    type( l2GPData_T ), intent(inout) :: l2gp
    integer, intent(in) :: l2FileHandle ! From swopen
    character (len=*), intent(in), optional :: swathName ! Defaults to l2gp%name

    ! Parameters

    character (len=*), parameter :: WR_ERR = 'Failed to write data field '

    ! Variables

    character (len=480) :: msr
    character (len=132) :: name     ! Either swathName or l2gp%name

    integer :: status
    integer :: start(3), stride(3), edge(3)
    integer :: swid

!  How to deal with status and columnTypes? swrfld fails
!  with char data on Linux
!  Have recourse to ints2Strings and strings2Ints if USEINTS4STRINGS
!    character (LEN=8), allocatable :: the_status_buffer(:)
!    character (LEN=L2GPNameLen), allocatable :: the_status_buffer(:)
    integer, allocatable, dimension(:,:) :: string_buffer

    if ( present(swathName) ) then
       name=swathName
    else
       name=l2gp%name
    end if
    ! Write data to the fields

    start = 0
    stride = 1
    edge(1) = l2gp%nFreqs
    edge(2) = l2gp%nLevels
    edge(3) = l2gp%nTimes
    swid = swattach (l2FileHandle, name)
    if ( l2gp%nFreqs > 0 ) then
       ! Value and Precision are 3-D fields

       status = swwrfld(swid, DATA_FIELD1, start, stride, edge, &
            & RESHAPE(l2gp%l2gpValue, (/SIZE(l2gp%l2gpValue)/)) )
       if ( status == -1 ) then
          msr = WR_ERR // DATA_FIELD1
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       status = swwrfld(swid, DATA_FIELD2, start, stride, edge, &
            & RESHAPE(real(l2gp%l2gpPrecision), (/SIZE(l2gp%l2gpPrecision)/)) )
       if ( status == -1 ) then
          msr = WR_ERR // DATA_FIELD2
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

    else if ( l2gp%nLevels > 0 ) then
       ! Value and Precision are 2-D fields

       status = swwrfld( swid, DATA_FIELD1, start(2:3), stride(2:3), &
            edge(2:3), real(l2gp%l2gpValue(1,:,:)) )

       if ( status == -1 ) then
          msr = WR_ERR // DATA_FIELD1
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       status = swwrfld( swid, DATA_FIELD2, start(2:3), stride(2:3), &
            edge(2:3), real(l2gp%l2gpPrecision(1,:,:) ))
       if ( status == -1 ) then
          msr = WR_ERR // DATA_FIELD2
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    else

       ! Value and Precision are 1-D fields
       status = swwrfld( swid, DATA_FIELD1, start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpValue(1,1,:) ))
       if ( status == -1 ) then
          msr = WR_ERR // DATA_FIELD1
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       status = swwrfld( swid, DATA_FIELD2, start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpPrecision(1,1,:) ))
       if ( status == -1 ) then
          msr = WR_ERR // DATA_FIELD2
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    end if

    ! 1-D status & quality fields

   if(USEINTS4STRINGS) then
      allocate(string_buffer(1,l2gp%nTimes))
      call strings2Ints(l2gp%status, string_buffer)
      status = swwrfld(swid, DATA_FIELD3, start(3:3), stride(3:3), edge(3:3), &
           string_buffer)
      if ( status == -1 ) then
         msr = WR_ERR // DATA_FIELD3
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
      deallocate(string_buffer)
   else
      status = swwrfld(swid, DATA_FIELD3, start(3:3), stride(3:3), edge(3:3), &
           l2gp%status)
      if ( status == -1 ) then
         msr = WR_ERR // DATA_FIELD3
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
   end if
    !  l2gp%quality = 0 !??????? Why was this here !??? NJL
    status = swwrfld(swid, DATA_FIELD4, start(3:3), stride(3:3), edge(3:3), &
         real(l2gp%quality))
    if ( status == -1 ) then
       msr = WR_ERR // DATA_FIELD4
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

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

    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    chunktimes=1
    chunkfreqs=1 ! better as nFreqs, but I have yet to see a case with nfreqs>1
    if(present(nLevels))then
       chunklevels=nLevels
    else
       chunklevels=5
    endif
    
    ! Create the swath within the file
    !print*,"Creating swath called ",name
    swid = HE5_SWcreate(L2FileHandle, trim(name))
    !print*,"Swath ",name,"has SW id :",swid
    if ( swid == -1 ) then
       msr = 'Failed to create swath ' // trim(name)
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    ! Define dimensions

    ! Defining special "unlimited dimension called UNLIM
    !print*,"Defined Unlim with size", H5S_UNLIMITED
    status = HE5_SWdefdim(swid, UNLIM, H5S_UNLIMITED)

    !print*,"Defining dimension ", DIM_NAME1," with size",l2gp%nTimes
    status = HE5_SWdefdim(swid, DIM_NAME1, l2gp%nTimes)
    if ( status == -1 ) then
       msr = DIM_ERR // DIM_NAME1
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    if ( l2gp%nLevels > 0 ) then
      !print*,"Defining dimension ", DIM_NAME2," with size",l2gp%nLevels
       status = HE5_SWdefdim(swid, DIM_NAME2, l2gp%nLevels)
       if ( status == -1 ) then
          msr = DIM_ERR // DIM_NAME2
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    end if

    if ( l2gp%nFreqs > 0 ) then
       print*,"Defining dimension ", DIM_NAME3," with size",l2gp%nFreqs
       status = HE5_SWdefdim(swid, DIM_NAME3, l2gp%nFreqs)
       if ( status == -1 ) then
          msr = DIM_ERR // DIM_NAME3
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    end if

    ! Define horizontal geolocation fields using above dimensions

!    print*,"Defining geolocation field ",GEO_FIELD1," of dim. ", DIM_NAME1
!    print*,"... and of type ",H5T_NATIVE_FLOAT
    chunk_rank=1
    chunk_dims=1
    chunk_dims(1)=CHUNKTIMES
    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
!    print*,"Set chunking -- status=",status
    status = HE5_SWdefgfld(swid, GEO_FIELD1, DIM_NAME1,MAX_DIML1,&
         H5T_NATIVE_FLOAT , 0)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD1
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if
!    print*,"Defined geolocation field ",GEO_FIELD1,"of dim.", DIM_NAME1
!    print*,"... and of type ",H5T_NATIVE_FLOAT

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)

    status = HE5_SWdefgfld(swid, GEO_FIELD2, DIM_NAME1, MAX_DIML1,&
         H5T_NATIVE_FLOAT, HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD2
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    status = HE5_SWdefgfld(swid, GEO_FIELD3, DIM_NAME1, MAX_DIML1, &
    H5T_NATIVE_DOUBLE, HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD3
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    status = HE5_SWdefgfld(swid, GEO_FIELD4, DIM_NAME1,MAX_DIML1,&
         H5T_NATIVE_FLOAT, HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD4
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    status = HE5_SWdefgfld(swid, GEO_FIELD5, DIM_NAME1, MAX_DIML1, &
         H5T_NATIVE_FLOAT, HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD5
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    status = HE5_SWdefgfld(swid, GEO_FIELD6, DIM_NAME1,MAX_DIML1,&
         H5T_NATIVE_FLOAT, HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD6
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    status = HE5_SWdefgfld(swid, GEO_FIELD7, DIM_NAME1, MAX_DIML1,&
    H5T_NATIVE_FLOAT,   HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD7
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    status = HE5_SWdefgfld(swid, GEO_FIELD8, DIM_NAME1, MAX_DIML1,&
         H5T_NATIVE_FLOAT, HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = GEO_ERR // GEO_FIELD8
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    if ( l2gp%nLevels > 0 ) then

       status = HE5_SWdefgfld(swid, GEO_FIELD9, DIM_NAME2,MAX_DIML,&
            H5T_NATIVE_FLOAT, HDFE_NOMERGE)
       if ( status == -1 ) then
          msr = GEO_ERR // GEO_FIELD9
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    end if

    if ( l2gp%nFreqs > 0 ) then

       status = HE5_SWdefgfld(swid, GEO_FIELD10, DIM_NAME3,MAX_DIML,&
            H5T_NATIVE_FLOAT, HDFE_NOMERGE)
       if ( status == -1 ) then
          msr = GEO_ERR // GEO_FIELD10
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    end if

    ! Define data fields using above dimensions

    if ( (l2gp%nFreqs > 0) .and. (l2gp%nLevels > 0) ) then
       chunk_rank=3
       chunk_dims(1:3)=(/ CHUNKFREQS,CHUNKLEVELS,CHUNKTIMES /)
       status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)

       status = HE5_SWdefdfld(swid, DATA_FIELD1, DIM_NAME123, MAX_DIML123,&
       H5T_NATIVE_FLOAT,HDFE_NOMERGE)

       if ( status == -1 ) then
          msr = DAT_ERR // DATA_FIELD1 // ' for 3D quantity.'
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if


       status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
       status = HE5_SWdefdfld(swid, DATA_FIELD2, DIM_NAME123, MAX_DIML123,&
            H5T_NATIVE_FLOAT, HDFE_NOMERGE)

       if ( status == -1 ) then
          msr = DAT_ERR // DATA_FIELD2 // ' for 3D quantity.'
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if


    else if ( l2gp%nLevels > 0 ) then
       chunk_rank=2
       chunk_dims(1:7)=(/ CHUNKLEVELS,CHUNKTIMES,37,38,39,47,49/)
       status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
       print*,"Set chunking with status=",status
       print*,"chunking=",chunk_dims
       print*,"About to define 2-D extendible field"

       print*,"Calling SWdefdfld with args ",swid, DATA_FIELD1, &
            DIM_NAME12, MAX_DIML12, H5T_NATIVE_FLOAT, HDFE_NOMERGE
       status = HE5_SWdefdfld(swid, DATA_FIELD1, DIM_NAME12, MAX_DIML12, &
            H5T_NATIVE_FLOAT, HDFE_NOMERGE)
       print*,"Defined 2-D extendible field"
    
       if ( status == -1 ) then
          msr = DAT_ERR // DATA_FIELD1 //  ' for 2D quantity.'
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

           status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
           status = HE5_SWdefdfld(swid, DATA_FIELD2, DIM_NAME12, MAX_DIML12,&
            H5T_NATIVE_FLOAT,HDFE_NOMERGE)

       if ( status == -1 ) then
          msr = DAT_ERR // DATA_FIELD2 //  ' for 2D quantity.'
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

    else
       chunk_rank=1
       chunk_dims(1)=CHUNKTIMES
       status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    
       status = HE5_SWdefdfld(swid, DATA_FIELD1, DIM_NAME1,MAX_DIML1,&
            H5T_NATIVE_FLOAT,HDFE_NOMERGE)

       if ( status == -1 ) then
          msr = DAT_ERR // DATA_FIELD1 // ' for 1D quantity.'
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

       status = HE5_SWdefdfld(swid, DATA_FIELD2, DIM_NAME1, MAX_DIML1,&
       H5T_NATIVE_FLOAT, HDFE_NOMERGE)

       if ( status == -1 ) then
          msr = DAT_ERR // DATA_FIELD2 // ' for 1D quantity.'
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

    end if

!    print*,"Defining data field ",DATA_FIELD3,"of dim.", DIM_NAME1
!    print*,"... and of type ",H5T_NATIVE_CHAR

!    status = HE5_SWdefdfld(swid, DATA_FIELD3, DIM_NAME1,MAX_DIML1,&
!         H5T_NATIVE_CHAR, HDFE_NOMERGE)
!    IF ( status == -1 ) THEN
!       msr = DAT_ERR // DATA_FIELD3
!       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
!    END IF

!    print*,"Defined data field ",DATA_FIELD3,"of dim.", DIM_NAME1
    chunk_rank=1
    chunk_dims(1)=CHUNKTIMES
    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)

    status = HE5_SWdefdfld(swid, DATA_FIELD4, DIM_NAME1,MAX_DIML1,&
         H5T_NATIVE_FLOAT, HDFE_NOMERGE)
    if ( status == -1 ) then
       msr = DAT_ERR // DATA_FIELD4
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

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

    character (len=480) :: msr
    character (len=132) :: name ! Either swathName or l2gp%name
    
    integer :: status, swid,myOffset
    integer :: start(2), stride(2), edge(2)

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
    start = myOffset
    edge(1) = l2gp%nTimes
    !print*,"writeGeo Attached swath ",name," with SW ID=",swid
    !print*,"About to write latitude with offset=",myoffset
    status = HE5_SWwrfld(swid, GEO_FIELD1, start, stride, edge, &
         real(l2gp%latitude))
    !print*,"wrote latitude, maybe"
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD1
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = HE5_SWwrfld(swid, GEO_FIELD2, start, stride, edge, &
         real(l2gp%longitude))
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD2
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = HE5_SWwrfld(swid, GEO_FIELD3, start, stride, edge, &
         l2gp%time)
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD3
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if


    !print*,"writing ", REAL(l2gp%solarZenith)," as SZA"
    status = HE5_SWwrfld(swid, GEO_FIELD5, start, stride, edge, &
         real(l2gp%solarZenith))
    !print*,"just wrote ", REAL(l2gp%solarZenith)," as SZA"
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD5
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = HE5_SWwrfld(swid, GEO_FIELD4, start, stride, edge, &
        real(l2gp%solarTime))
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD4
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if


    status = HE5_SWwrfld(swid, GEO_FIELD6, start, stride, edge, &
         real(l2gp%losAngle))
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD6
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = HE5_SWwrfld(swid, GEO_FIELD7, start, stride, edge, &
         real(l2gp%geodAngle))
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD7
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    status = HE5_SWwrfld(swid, GEO_FIELD8, start, stride, edge, &
         l2gp%chunkNumber)
    if ( status == -1 ) then
       msr = WR_ERR // GEO_FIELD8
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    if ( l2gp%nLevels > 0 ) then
       edge(1) = l2gp%nLevels
       start(1)=0 ! needed because offset may have made this /=0
       status = HE5_SWwrfld(swid, GEO_FIELD9, start, stride, edge, &
            real(l2gp%pressures))
       if ( status == -1 ) then
          msr = WR_ERR // GEO_FIELD9
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    end if

    if ( l2gp%nFreqs > 0 ) then
       edge(1) = l2gp%nFreqs
       start(1)=0 ! needed because offset may have made this /=0
       l2gp%frequency = 0
       status = HE5_SWwrfld(swid, GEO_FIELD10, start, stride, edge, &
            real(l2gp%frequency))
       if ( status == -1 ) then
          msr = WR_ERR // GEO_FIELD10
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    end if

    ! Detach from the swath interface.  

    status = HE5_SWdetach(swid)
    !print*,"Detatched from swath -- error=",status
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach from swath interface' )
    end if

    !------------------------------------
  end subroutine OutputL2GP_writeGeo_hdf5
  !------------------------------------

  !----------------------------------------  OutputL2GP_writeData_hdf5  -----
  subroutine OutputL2GP_writeData_hdf5(l2gp, l2FileHandle, swathName,offset)

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

    !print*,"OutputL2GP_writeData -- name=",name
    ! Write data to the fields

    start = 0
    stride = 1
    start(3)=myOffset
    edge(1) = l2gp%nFreqs
    edge(2) = l2gp%nLevels
    edge(3) = l2gp%nTimes
    swid = HE5_SWattach (l2FileHandle, name)
    !print*," attached swath with swid=",swid," filehandle=",l2FileHandle
    if ( l2gp%nFreqs > 0 ) then
       !print*,"Writing 3D field"
       ! Value and Precision are 3-D fields
       status = HE5_SWwrfld(swid, DATA_FIELD1, start, stride, edge, &
            & reshape(l2gp%l2gpValue, (/size(l2gp%l2gpValue)/)) )
       if ( status == -1 ) then
          msr = WR_ERR // DATA_FIELD1
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       status = HE5_SWwrfld(swid, DATA_FIELD2, start, stride, edge, &
            & reshape(real(l2gp%l2gpPrecision), (/size(l2gp%l2gpPrecision)/)) )
       if ( status == -1 ) then
          msr = WR_ERR // DATA_FIELD2
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if

    else if ( l2gp%nLevels > 0 ) then
       !Print*,"Writing 2-d field"
       ! Value and Precision are 2-D fields
     ! print*,"About to write data with offset=",myOffset
      
       status = HE5_SWwrfld( swid, DATA_FIELD1, start(2:3), stride(2:3), &
            edge(2:3), real(l2gp%l2gpValue(1,:,:) ))
       !print*,"Status of write was ",status
       if ( status == -1 ) then
          msr = WR_ERR // DATA_FIELD1
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       status = HE5_SWwrfld( swid, DATA_FIELD2, start(2:3), stride(2:3), &
            edge(2:3), real(l2gp%l2gpPrecision(1,:,:) ))
       if ( status == -1 ) then
          msr = WR_ERR // DATA_FIELD2
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    else

       ! Value and Precision are 1-D fields
       !Print*,"Writing 1-D field"
       status = HE5_SWwrfld( swid, DATA_FIELD1, start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpValue(1,1,:) ))
       if ( status == -1 ) then
          msr = WR_ERR // DATA_FIELD1
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
       status = HE5_SWwrfld( swid, DATA_FIELD2, start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpPrecision(1,1,:) ))
       if ( status == -1 ) then
          msr = WR_ERR // DATA_FIELD2
          call MLSMessage ( MLSMSG_Error, ModuleName, msr )
       end if
    end if

    ! 1-D status & quality fields

    !HDF-EOS5 won't write a dataset of chars from FORTRAN

   if(USEINTS4STRINGS) then
      allocate(string_buffer(1,l2gp%nTimes))
      call strings2Ints(l2gp%status, string_buffer)
      status = HE5_swwrfld(swid,DATA_FIELD3,start(3:3),stride(3:3),edge(3:3),&
           string_buffer)
      if ( status == -1 ) then
         msr = WR_ERR // DATA_FIELD3
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
      deallocate(string_buffer)
    else
      !    status = HE5_SWwrfld(swid, DATA_FIELD3, start(3:3), stride(3:3),&
      !        edge(3:3), l2gp%status) ! 
      status=0
      print*,"Warning. Writing of status field disabled"
      if ( status == -1 ) then
         msr = WR_ERR // DATA_FIELD3
         call MLSMessage ( MLSMSG_Error, ModuleName, msr )
      end if
    end if
    !  l2gp%quality = 0 !??????? Why was this here !??? NJL
    !                   ! Beats me. Evil bug gnomes, probably. HCP 
    status = HE5_SWwrfld(swid, DATA_FIELD4, start(3:3), stride(3:3), edge(3:3), &
         real(l2gp%quality))
    if ( status == -1 ) then
       msr = WR_ERR // DATA_FIELD4
       call MLSMessage ( MLSMSG_Error, ModuleName, msr )
    end if

    !     Detach from the swath interface.

    status = HE5_SWdetach(swid)
    print*,"Detatched from swath -- error=",status
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach  from swath interface' )
    end if


    !-------------------------------------
  end subroutine OutputL2GP_writeData_hdf5
  !-------------------------------------

  ! --------------------------------------------------------------------------

  ! This subroutine is an amalgamation of the last three
  ! Should be renamed CreateAndWriteL2GPData
  subroutine WriteL2GPData(l2gp,l2FileHandle,swathName, hdfVersion)

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

    if (myhdfVersion == 4) then
      call OutputL2GP_createFile_hdf4 (l2gp, l2FileHandle, swathName)
      call OutputL2GP_writeGeo_hdf4 (l2gp, l2FileHandle, swathName)
      call OutputL2GP_writeData_hdf4 (l2gp, l2FileHandle, swathName)
    else
      call OutputL2GP_createFile_hdf5 (l2gp, l2FileHandle, swathName)
      call OutputL2GP_writeGeo_hdf5 (l2gp, l2FileHandle, swathName)
      call OutputL2GP_writeData_hdf5 (l2gp, l2FileHandle, swathName)
    endif

  end subroutine WriteL2GPData
  !-------------------------------------------------------------


  subroutine AppendL2GPData(l2gp,l2FileHandle,swathName,offset, hdfVersion)
    ! sticks l2gp into the swath swathName in the file pointed at by
    ! l2FileHandle,starting at the profile number "offset" (First profile
    ! in the file has offset==0). If this runs off the end of the swath, 
    ! it is lengthend automagically. 
    ! Arguments

    integer, intent(IN) :: l2FileHandle ! From swopen
    type (L2GPData_T), intent(INOUT) :: l2gp
    character (LEN=*), optional, intent(IN) ::swathName!default->l2gp%swathName
    integer,intent(IN),optional::offset
    integer, optional, intent(in) :: hdfVersion
    !----Local variable
    integer::myOffset
    ! Local
    integer :: myhdfVersion

    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L2GPDEFAULT_HDFVERSION
    endif

    ! ----Executable code---
    if (present(offset)) then
       myOffset=offset
    else
       myOffset=0
    endif
    if (myhdfVersion == 4) then
      call OutputL2GP_writeGeo_hdf4 (l2gp, l2FileHandle, swathName)
      call OutputL2GP_writeData_hdf4 (l2gp, l2FileHandle, swathName)
    else
      call OutputL2GP_writeGeo_hdf5 (l2gp, l2FileHandle, swathName,myOffset)
      call OutputL2GP_writeData_hdf5 (l2gp, l2FileHandle, swathName,myOffset)
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
    integer :: i, ierr
    logical :: myColumnsOnly
    integer :: MYDETAILS

    ! Executable code
    myDetails = 1
    if ( present(details) ) myDetails = details
    
    if( present(ColumnsOnly)) then
      myColumnsOnly = ColumnsOnly
    else
      myColumnsOnly = .false.
    endif

      if ( myColumnsOnly .and. l2gp%nLevels > 1 ) return

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
      call dump ( l2gp%l2gpValue, 'L2GPValue:' )
      
      call dump ( l2gp%l2gpPrecision, 'L2GPPrecision:' )
      
      !    call dump ( l2gp%status, 'Status:' )
      
      call dump ( l2gp%quality, 'Quality:' )
      
  end subroutine Dump_L2GP
    

  !=============================================================================
end module L2GPData
!=============================================================================

!
! $Log$
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
