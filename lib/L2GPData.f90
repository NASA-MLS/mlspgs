! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module L2GPData                 ! Creation, manipulation and I/O for L2GP Data
  !=============================================================================

  use Allocate_Deallocate, only: Allocate_test, Deallocate_test
  use DUMP_0, only: DUMP
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
       & MLSMSG_Error, MLSMSG_Warning
  use MLSCommon, only: R8
  use MLSStrings, only: ints2Strings, strings2Ints
  use Hdf, only: DFNT_CHAR8, DFNT_FLOAT32, DFNT_INT32, DFNT_FLOAT64
  use HDFEOS!, only: SWATTACH, SWCREATE, SWDEFDFLD, SWDEFDIM, SWDEFGFLD, &
     !& SWDETACH
  use OUTPUT_M, only: OUTPUT
  use SWAPI, only: SWWRFLD, SWRDFLD
  use STRING_TABLE, only: DISPLAY_STRING
 
  implicit none

  private
  public :: L2GPData_T
  public :: L2GPNameLen
  public :: AddL2GPToDatabase,  DestroyL2GPContents,  DestroyL2GPDatabase, &
    & Dump, Dump_L2GP,  ExpandL2GPDataInPlace,  OutputL2GP_createFile, &
    & OutputL2GP_writeData,  OutputL2GP_writeGeo,  ReadL2GPData, &
    & SetupNewL2GPRecord,  WriteL2GPData
!  INTEGER :: SWRDFLD
!  EXTERNAL SWRDFLD !Should USE SWAPI
  !---------------------------- RCS Ident Info -------------------------------
  character(len=*), private, parameter :: IdParm = &
    & "$Id$"
  character(len=len(idparm)), private :: Id = idParm
  character(len=*), private, parameter :: ModuleName = &
       & "$RCSfile$"
  !---------------------------------------------------------------------------

  interface DUMP
    module procedure DUMP_L2GP
  end interface


  ! This module defines datatypes and gives basic routines for storing and
  ! manipulating L2GP data.

  ! First some local parameters

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

  character (len=*), parameter :: DIM_NAME12 = 'nLevels,nTimes' ! In Fortran order?!!
  character (len=*), parameter :: DIM_NAME123 = 'nFreqs,nLevels,nTimes' ! as above

  integer, parameter :: HDFE_AUTOMERGE = 1     ! MERGE FIELDS WITH SHARE DIM
  integer, parameter :: HDFE_NOMERGE = 0       ! don't merge
  ! The following, if true, says to encode strings as ints
  ! before swapi write; also decode ints to strings after read
  ! otherwise, try swapi read, write directly with strings
  logical, parameter :: USEINTS4STRINGS = .false.  
  
  ! This datatype is the main one, it simply defines one l2gp swath
  ! It is used for 
  ! (1) normal swaths, which have geolocations along nTimes
  !     vertical coordinates called "surfaces" dimensioned by pressures
  !     and at nTimes horizontal "instances" dimensioned by 
  !     latitudes and longitudes
  ! (2) column abundances which are  integrated amounts and therefore 
  !     the vertical coordinate is suppressed; instead we set nLevels=1
  !     and store the tropopause in pressures(1)

  typE L2GPData_T

     ! First some variables we use for book keeping

     character (len=l2gpnamelen) :: name ! Typically the swath name.
     integer :: nameIndex       ! Used by the parser to keep track of the data

     ! Now the dimensions of the data

     integer :: nTimes          ! Total number of profiles
     integer :: nLevels         ! Total number of surfaces (=1 for column abundances)
     integer :: nFreqs          ! Number of frequencies in breakdown

     ! Now we store the geolocation fields, first the vertical one:
     ! (The following has the tropopause if the swath is a column abundance)
     real (r8), pointer, dimension(:) :: pressures => NULL() ! Vertical coords (nLevels)

     ! Now the horizontal geolocation information. Dimensioned (nTimes)
     real (r8), pointer, dimension(:) :: latitude => NULL()
     real (r8), pointer, dimension(:) :: longitude => NULL()
     real (r8), pointer, dimension(:) :: solarTime => NULL()
     real (r8), pointer, dimension(:) :: solarZenith => NULL()
     real (r8), pointer, dimension(:) :: losAngle => NULL()
     real (r8), pointer, dimension(:) :: geodAngle => NULL()
     real (r8), pointer, dimension(:) :: time => NULL()

     integer, pointer, dimension(:) :: chunkNumber => NULL()

     ! Now we store the `frequency' geolocation field

     real (r8), pointer, dimension(:) :: frequency => NULL()
     !        dimensioned (nFreqs)

     ! Finally we store the data fields

     real (r8), pointer, dimension(:,:,:) :: l2gpValue => NULL()
     real (r8), pointer, dimension(:,:,:) :: l2gpPrecision => NULL()
     ! dimensioned (nFreqs, nLevels, nTimes)

     character (len=1), pointer, dimension(:) :: status => NULL()
     real (r8), pointer, dimension(:) :: quality => NULL()
     ! Both the above dimensioned (nTimes)

  end type L2GPData_T

contains ! =====     Public Procedures     =============================

  !------------------------------------------  SetupNewL2GPRecord  -----
  subroutine SetupNewL2GPRecord ( l2gp, nFreqs, nLevels, nTimes)

    ! This routine sets up the arrays for an l2gp datatype.

    ! Dummy arguments
    type (L2GPData_T), intent(inout)  :: l2gp
    integer, intent(in), optional :: nFreqs            ! Dimensions
    integer, intent(in), optional :: nLevels           ! Dimensions
    integer, intent(in), optional :: nTimes            ! Dimensions

    ! Local variables
    integer :: freqsArrayLen, status, surfsArrayLen
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
    call SetupNewL2GPRecord( l2gp, nFreqs=l2gp%nFreqs, nLevels=l2gp%nLevels, &
      & nTimes=myNTimes )

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

  ! -------------------------------------------------------------------------

  subroutine ReadL2GPData(L2FileHandle, swathname, l2gp, numProfs, &
       firstProf, lastProf)
    !------------------------------------------------------------------------

    ! This routine reads an L2GP file, returning a filled data structure and the !
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
  end subroutine ReadL2GPData
  !-----------------------------

  ! --------------------------------------  OutputL2GP_createFile  -----
  subroutine OutputL2GP_createFile (l2gp, L2FileHandle, swathName)

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
  end subroutine OutputL2GP_createFile
  !--------------------------------------

  !-----------------------------------------  OutputL2GP_writeGeo  -----
  subroutine OutputL2GP_writeGeo (l2gp, l2FileHandle, swathName)

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
  end subroutine OutputL2GP_writeGeo
  !------------------------------------

  !----------------------------------------  OutputL2GP_writeData  -----
  subroutine OutputL2GP_writeData(l2gp, l2FileHandle, swathName)

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
  end subroutine OutputL2GP_writeData
  !-------------------------------------

  ! --------------------------------------------------------------------------

  ! This subroutine is an amalgamation of the last three

  subroutine WriteL2GPData(l2gp,l2FileHandle,swathName)

    ! Arguments

    integer, intent(in) :: l2FileHandle ! From swopen
    type (l2GPData_T), intent(inout) :: l2gp
    character (len=*), optional, intent(in) :: swathName ! (defaults to l2gp%swathName)

    ! Exectuable code

    call OutputL2GP_createFile (l2gp, l2FileHandle, swathName)
    call OutputL2GP_writeGeo (l2gp, l2FileHandle, swathName)
    call OutputL2GP_writeData (l2gp, l2FileHandle, swathName)

  end subroutine WriteL2GPData

  ! ------------------------------------------ DUMP_L2GP ------------

  subroutine Dump_L2GP ( L2gp, Name, ColumnsOnly )

    ! Dummy arguments
    type (l2gpData_T), intent(in) ::          L2GP(:)
    character(len=*), intent(in), optional :: Name
    logical, intent(in), optional ::          ColumnsOnly ! if true, dump only with columns

    ! Local variables
    integer :: i
    logical :: myColumnsOnly
    
    if( present(ColumnsOnly)) then
      myColumnsOnly = ColumnsOnly
    else
      myColumnsOnly = .false.
    endif

    if ( present(name) ) call output ( name, advance='yes' )
    do i = 1, size(l2gp)
      call output ( 'L2GP Data: ')
      if(l2gp(i)%nameIndex > 0) then
        call display_string ( l2gp(i)%nameIndex, advance='yes' )
      else
        call output ( '(the nameIndex is 0) ', advance='yes')
      endif
      call output ( 'nTimes: ')
      call output ( l2gp(i)%nTimes, 5)
      call output ( '  nLevels: ')
      call output ( l2gp(i)%nLevels, 3)
      call output ( '  nFreqs: ')
      call output ( l2gp(i)%nFreqs, 3, advance='yes')
      if(myColumnsOnly .and. l2gp(i)%nLevels>1) CYCLE
      call dump ( l2gp(i)%pressures, 'Pressures:' )
      
      call dump ( l2gp(i)%latitude, 'Latitude:' )
      
      call dump ( l2gp(i)%longitude, 'Longitude:' )
      
      call dump ( l2gp(i)%solarTime, 'SolarTime:' )
      
      call dump ( l2gp(i)%solarZenith, 'SolarZenith:' )
      
      call dump ( l2gp(i)%losAngle, 'LOSAngle:' )
      
      call dump ( l2gp(i)%geodAngle, 'geodAngle:' )
      
      call dump ( l2gp(i)%time, 'Time:' )
      
      call dump ( l2gp(i)%chunkNumber, 'ChunkNumber:' )
      
      if ( associated(l2gp(i)%frequency) ) &
        & call dump ( l2gp(i)%frequency, 'Frequencies:' )
      
      call dump ( l2gp(i)%l2gpValue, 'L2GPValue:' )
      
      call dump ( l2gp(i)%l2gpPrecision, 'L2GPPrecision:' )
      
      !    call dump ( l2gp(i)%status, 'Status:' )
      
      call dump ( l2gp(i)%quality, 'Quality:' )
      
    end do
  end subroutine Dump_L2GP
    
  !=============================================================================
end module L2GPData
!=============================================================================

!
! $Log$
! Revision 2.36  2001/08/16 23:17:42  vsnyder
! Remove continuation ampersands within two lines
!
! Revision 2.35  2001/08/06 18:39:03  pwagner
! Now dumps column-related components, too
!
! Revision 2.34  2001/08/03 23:13:52  pwagner
! Began testing; at least now exits normally again
!
! Revision 2.33  2001/08/03 00:02:26  pwagner
! Ncolumns now a component of data type; swapi via ints not strings
!
! Revision 2.32  2001/08/02 00:17:56  pwagner
! Added column components; untested
!
! Revision 2.31  2001/06/13 20:36:15  vsnyder
! Commented out reading l2gp%status.  It appears that reading characters from
! HDF can't be made to work, due to fundamental design flaws in HDF.
!
! Revision 2.30  2001/06/06 17:28:13  pwagner
! the_status_buffer allows reading l2gp%status
!
! Revision 2.29  2001/05/03 23:59:56  vsnyder
! Trying to find out why realProf is undefined
!
! Revision 2.28  2001/04/24 16:28:41  livesey
! Cosmetic changes only, except comment out dubious l2gp%quality = 0 statement.
!
! Revision 2.27  2001/04/20 03:00:59  livesey
! Made L2GPNameLen public
!
! Revision 2.26  2001/04/20 02:05:09  vsnyder
! Cosmetic changes: Default visibility is now private
!
! Revision 2.25  2001/03/20 01:44:25  livesey
! Fixed bug, was outputting chunkNumber as real!
!
! Revision 2.24  2001/03/12 20:25:04  vsnyder
! Improve dump_l2gp.  Nullify components of l2gp before calling
! SetupNewL2GPRecord from ExpandL2GPDataInPlace, so as not to clobber the
! component pointers carefully copied to tempL2gp.
!
! Revision 2.23  2001/03/01 18:37:51  livesey
! Added dumper routine
!
! Revision 2.22  2001/02/22 21:54:22  livesey
! Added initialisation to NULL() for pointer components of L2GPData_T
!
! Revision 2.21  2001/02/15 18:23:20  livesey
! Got it right this time!
!
! Revision 2.20  2001/02/15 18:20:15  livesey
! Had to comment out reading of status in ReadL2GPData to make it work, logged PR.
!
! Revision 2.19  2001/02/14 23:39:40  livesey
! Made numProfs argument optional(intent out) for ReadL2GPData
!
! Revision 2.18  2001/02/13 00:55:35  livesey
! Removed another print statement!
!
! Revision 2.17  2001/02/13 00:51:31  livesey
! Fixed bug, copy pressure information in ExpandL2GPDataInPlace
!
! Revision 2.16  2001/02/09 18:38:04  livesey
! Even more print statemets removed!
!
! Revision 2.15  2001/02/09 17:51:01  livesey
! Removed some print statements.
!
! Revision 2.14  2001/02/09 17:45:15  livesey
! Another fix to dimension ordering.
!
! Revision 2.13  2001/02/08 01:06:08  livesey
! Bug fix in ExpandL2GPDataInPlace
!
! Revision 2.12  2001/02/05 23:58:22  pwagner
! Uses swrdfld from swapi
!
! Revision 2.11  2001/02/03 00:04:26  pwagner
! Uncommented EXTERNAL SWRDFLD until swapi done
!
! Revision 2.10  2001/02/02 17:15:02  livesey
! Changed order of DIM_NAME12 and DIM_NAME123, this should be right now!
!
! Revision 2.9  2001/01/29 18:17:49  livesey
! Changed status, quality and frequency to upper case first character.
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
