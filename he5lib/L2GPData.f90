! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L2GPData                 ! Creation, manipulation and I/O for L2GP Data
!=============================================================================
  USE Allocate_Deallocate, ONLY: Allocate_test, Deallocate_test
  USE DUMP_0, only: DUMP
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
       & MLSMSG_Error, MLSMSG_Warning
  USE MLSCommon, ONLY: R8

  use HDF5_params
  USE HDFEOS5
  USE HE5_SWAPI 

  USE OUTPUT_M, only: OUTPUT ! Added as HDF4 version uses it
  USE STRING_TABLE, only: DISPLAY_STRING

  IMPLICIT NONE
  !INTEGER:: HE5_SWRDFLD
  !EXTERNAL HE5_SWRDFLD !Should USE HE5_SWAPI
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER(len=256), PRIVATE :: Id = &
       & "$Id$"
  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName = &
       & "$RCSfile$"
  !---------------------------------------------------------------------------


  interface DUMP !And this does WTF? 
    module procedure DUMP_L2GP
  end interface

  ! This module defines datatypes and gives basic routines for storing and
  ! manipulating L2GP data.

  ! First some local parameters

  INTEGER, PARAMETER :: L2GPNameLen = 80

   CHARACTER (len=*), PARAMETER :: DATA_FIELD1 = 'L2gpValue'
   CHARACTER (len=*), PARAMETER :: DATA_FIELD2 = 'L2gpPrecision'
   CHARACTER (len=*), PARAMETER :: DATA_FIELD3 = 'Status'
   CHARACTER (len=*), PARAMETER :: DATA_FIELD4 = 'Quality'

   CHARACTER (len=*), PARAMETER :: GEO_FIELD1 = 'Latitude'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD2 = 'Longitude'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD3 = 'Time'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD4 = 'LocalSolarTime'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD5 = 'SolarZenithAngle'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD6 = 'LineOfSightAngle'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD7 = 'OrbitGeodeticAngle'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD8 = 'ChunkNumber'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD9 = 'Pressure'
   CHARACTER (len=*), PARAMETER :: GEO_FIELD10= 'frequency'

   CHARACTER (len=*), PARAMETER :: DIM_NAME1 = 'nTimes'
   CHARACTER (len=*), PARAMETER :: DIM_NAME2 = 'nLevels'
   CHARACTER (len=*), PARAMETER :: DIM_NAME3 = 'nFreqs'
   CHARACTER (len=*), PARAMETER :: DIM_NAME12 = 'nLevels,nTimes'
   CHARACTER (len=*), PARAMETER :: DIM_NAME123 = 'nFreqs,nLevels,nTimes'
   ! These are for the new max_dimlist parameter added to  SWdefgfld.
   ! this one is for non-extendible dimensions
   CHARACTER (len=*), PARAMETER :: MAX_DIML = ' '
   CHARACTER (len=*), PARAMETER :: UNLIM = 'Unlim'
   ! This is for cases where the time dimension is extendible
   CHARACTER (len=*), PARAMETER :: MAX_DIML1 = UNLIM
   CHARACTER (len=*), PARAMETER :: MAX_DIML12 = 'nLevels,Unlim'
   CHARACTER (len=*), PARAMETER :: MAX_DIML123 = 'nFreqs,nLevels,Unlim'

!   INTEGER,PARAMETER::CHUNKFREQS=13,CHUNKLEVELS=17,CHUNKTIMES=9,CHUNK4=1

   INTEGER, PARAMETER :: HDFE_AUTOMERGE = 1     ! MERGE FIELDS WITH SHARE DIM
   INTEGER, PARAMETER :: HDFE_NOMERGE = 0       ! don't merge

  ! This datatype is the main one, it simply defines one l2gp swath

  TYPE L2GPData_T

     ! First some variables we use for book keeping

     CHARACTER (LEN=L2GPNameLen) :: name ! Typically the swath name.
     INTEGER :: nameIndex       ! Used by the parser to keep track of the data

     ! Now the dimensions of the data

     INTEGER :: nTimes          ! Total number of profiles
     INTEGER :: nLevels         ! Total number of surfaces
     INTEGER :: nFreqs          ! Number of frequencies in breakdown

     ! Now we store the geolocation fields, first the vertical one:
     REAL (r8), POINTER, DIMENSION(:) :: pressures=>NULL() ! Vertical coords (nLevels)

     ! Now the horizontal geolocation information. Dimensioned (nTimes)
     REAL (r8), POINTER, DIMENSION(:) :: latitude => NULL()
     REAL (r8), POINTER, DIMENSION(:) :: longitude => NULL()
     REAL (r8), POINTER, DIMENSION(:) :: solarTime => NULL()
     REAL (r8), POINTER, DIMENSION(:) :: solarZenith => NULL()
     REAL (r8), POINTER, DIMENSION(:) :: losAngle => NULL()
     REAL (r8), POINTER, DIMENSION(:) :: geodAngle => NULL()
     REAL (r8), POINTER, DIMENSION(:) :: time => NULL()

     INTEGER, POINTER, DIMENSION(:) :: chunkNumber=>NULL()

     ! Now we store the `frequency' geolocation field

     REAL (r8), POINTER, DIMENSION(:) :: frequency=>NULL()
     !        dimensioned (nFreqs)

     ! Finally we store the data fields

     REAL (r8), POINTER, DIMENSION(:,:,:) :: l2gpValue=>NULL()
     REAL (r8), POINTER, DIMENSION(:,:,:) :: l2gpPrecision=>NULL()
     ! dimensioned (nFreqs, nLevels, nTimes)

     CHARACTER (len=1), POINTER, DIMENSION(:) :: status=>NULL()
     !                (status is a reserved word in F90)
     REAL (r8), POINTER, DIMENSION(:) :: quality=>NULL()
     ! Both the above dimensioned (nTimes)

  END TYPE L2GPData_T

CONTAINS ! =====     Public Procedures     =============================

  !------------------------------------------  SetupNewL2GPRecord  -----
  SUBROUTINE SetupNewL2GPRecord ( l2gp, nFreqs, nLevels, nTimes)

    ! This routine sets up the arrays for an l2gp datatype.

    ! Dummy arguments
    TYPE (L2GPData_T), INTENT(out)  :: l2gp
    INTEGER, INTENT(in), OPTIONAL :: nFreqs, nLevels, nTimes ! Dimensions

    ! Local variables
    INTEGER :: freqsArrayLen, status, surfsArrayLen
    INTEGER :: useNFreqs, useNLevels, useNTimes

    IF (PRESENT(nFreqs)) THEN
       useNFreqs=nFreqs
    ELSE
       useNFreqs=0
    ENDIF

    IF (PRESENT(nLevels)) THEN
       useNLevels=nLevels
    ELSE
       useNLevels=0
    ENDIF

    IF (PRESENT(nTimes)) THEN
       useNTimes=nTimes
    ELSE
       useNTimes=0              ! Default to empty l2gp
    ENDIF

    ! Store the dimensionality

    l2gp%nTimes = useNTimes
    l2gp%nLevels = useNLevels
    l2gp%nFreqs = useNFreqs

    ! But allocate to at least one for times, freqs

    useNTimes=MAX(useNTimes,1)
    useNLevels=MAX(useNLevels,1)
    useNFreqs=MAX(useNFreqs,1)    

    ! Allocate the frequency coordinate

    CALL allocate_test ( l2gp%pressures, useNLevels, "l2gp%pressures", &
         & ModuleName )

    ! Allocate the vertical coordinate

    CALL allocate_test ( l2gp%frequency, useNFreqs, "l2gp%frequency", ModuleName)

    ! Allocate the horizontal coordinates

    CALL allocate_test(l2gp%latitude,   useNTimes, "l2gp%latitude",   ModuleName)
    CALL allocate_test(l2gp%longitude,  useNTimes, "l2gp%longitude",  ModuleName)
    CALL allocate_test(l2gp%solarTime,  useNTimes, "l2gp%solarTime",  ModuleName)
    CALL allocate_test(l2gp%solarZenith,useNTimes, "l2gp%solarZenith",ModuleName)
    CALL allocate_test(l2gp%losAngle,   useNTimes, "l2gp%losAngle",   ModuleName)
    CALL allocate_test(l2gp%geodAngle,  useNTimes, "l2gp%geodAngle",  ModuleName)
    CALL allocate_test(l2gp%time,       useNTimes, "l2gp%time",       ModuleName)
    CALL allocate_test(l2gp%chunkNumber,useNTimes, "l2gp%chunkNumber",ModuleName)

    ! Allocate the data fields

    CALL allocate_test(l2gp%l2gpValue,useNFreqs,useNLevels,useNTimes,"l2gp%l2gpValue", ModuleName)
    CALL allocate_test(l2gp%l2gpPrecision,useNFreqs,useNLevels,useNTimes,"l2gp%l2gpPrecision", ModuleName)

    CALL allocate_test(l2gp%status, useNTimes,"l2gp%status", ModuleName)
    CALL allocate_test(l2gp%quality,useNTimes,"l2gp%quality",ModuleName)

  END SUBROUTINE SetupNewL2GPRecord

  !-----------------------------------------  DestroyL2GPContents  -----
  SUBROUTINE DestroyL2GPContents ( L2GP )

    ! This routine deallocates all the arrays allocated above.

    ! Dummy arguments
    TYPE (L2GPData_T), INTENT(inout) :: L2GP
    ! Local variables

    INTEGER status

    ! Executable code

    CALL deallocate_test ( l2gp%pressures,    "l2gp%pressures",    ModuleName )
    CALL deallocate_test ( l2gp%latitude,     "l2gp%latitude",     ModuleName )
    CALL deallocate_test ( l2gp%longitude,    "l2gp%longitude",    ModuleName )
    CALL deallocate_test ( l2gp%solarTime,    "l2gp%solarTime",    ModuleName )
    CALL deallocate_test ( l2gp%solarZenith,  "l2gp%solarZenith",  ModuleName )
    CALL deallocate_test ( l2gp%losAngle,     "l2gp%losAngle",     ModuleName )
    CALL deallocate_test ( l2gp%losAngle,     "l2gp%losAngle",     ModuleName )
    CALL deallocate_test ( l2gp%geodAngle,    "l2gp%geodAngle",    ModuleName )
    CALL deallocate_test ( l2gp%chunkNumber,  "l2gp%chunkNumber",  ModuleName )
    CALL deallocate_test ( l2gp%time,         "l2gp%time",         ModuleName )
    CALL deallocate_test ( l2gp%frequency,    "l2gp%frequency",    ModuleName )
    CALL deallocate_test ( l2gp%l2gpValue,    "l2gp%l2gpValue",    ModuleName )
    CALL deallocate_test ( l2gp%l2gpPrecision,"l2gp%l2gpPrecision",ModuleName )
    CALL deallocate_test ( l2gp%status,       "l2gp%status",       ModuleName )
    CALL deallocate_test ( l2gp%quality,      "l2gp%quality",      ModuleName )
    l2gp%nTimes = 0
    l2gp%nLevels = 0
    l2gp%nFreqs = 0
  END SUBROUTINE DestroyL2GPContents

  !---------------------------------------  ExpandL2GPDataInPlace  -----
  SUBROUTINE ExpandL2GPDataInPlace ( l2gp, newNTimes )

    ! This subroutine expands an L2GPData_T in place allowing the user to add
    ! more profiles to it.

    ! Dummy arguments
    TYPE (L2GPData_T), INTENT(inout) :: l2gp
    INTEGER, INTENT(in) :: newNTimes

    ! Local variables
    INTEGER :: status                   ! From ALLOCATE
    TYPE (L2GPData_T) :: tempL2gp       ! For copying data around

    ! Executable code
    ! First do a sanity check

    IF ( newNTimes<l2gp%nTimes ) &
         & CALL MLSMessage ( MLSMSG_Error, ModuleName, &
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
    CALL SetupNewL2GPRecord( l2gp, nFreqs=l2gp%nFreqs, nLevels=l2gp%nLevels, &
      & nTimes=newNTimes)
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
    CALL DestroyL2GPContents(templ2gp)

  END SUBROUTINE ExpandL2GPDataInPlace

  !-------------------------------------------  AddL2GPToDatabase  -----
  INTEGER FUNCTION AddL2GPToDatabase( DATABASE, ITEM )

    ! This function adds an l2gp data type to a database of said types,
    ! creating a new database if it doesn't exist.  The result value is
    ! the size -- where L2gp is put.

    ! Dummy arguments
    TYPE (l2gpdata_t), DIMENSION(:), POINTER :: DATABASE
    TYPE (l2gpdata_t), INTENT(in) :: ITEM

    ! Local variables
    TYPE (L2GPData_T), DIMENSION(:), POINTER :: tempDatabase
    !This include causes real trouble if you are compiling in a different 
    !directory.
    INCLUDE "addItemToDatabase.f9h" 

    AddL2GPToDatabase = newSize
  END FUNCTION AddL2GPToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template database

  SUBROUTINE DestroyL2GPDatabase ( DATABASE )

    ! Dummy argument
    TYPE (L2GPData_T), DIMENSION(:), POINTER :: DATABASE

    ! Local variables
    INTEGER :: l2gpIndex, status

    IF ( ASSOCIATED(database)) THEN
       DO l2gpIndex = 1, SIZE(database)
          CALL DestroyL2GPContents ( database(l2gpIndex) )
       END DO
       DEALLOCATE ( database, stat=status )
       IF ( status /= 0 ) CALL MLSMessage ( MLSMSG_Error, ModuleName, &
            & MLSMSG_deallocate // "database" )
    END IF
  END SUBROUTINE DestroyL2GPDatabase

  ! -------------------------------------------------------------------------

  SUBROUTINE ReadL2GPData(L2FileHandle, swathname, l2gp, numProfs, &
       firstProf, lastProf)
    !------------------------------------------------------------------------

    ! This routine reads an L2GP file, returning a filled data structure and the !
    ! number of profiles read.

    ! Arguments

    CHARACTER (LEN=*), INTENT(IN) :: swathname ! Name of swath
    INTEGER, INTENT(IN) :: L2FileHandle ! Returned by swopen
    INTEGER, INTENT(IN), OPTIONAL :: firstProf, lastProf ! Defaults to first and last
    TYPE( L2GPData_T ), INTENT(OUT) :: l2gp ! Result
    INTEGER, INTENT(OUT),OPTIONAL :: numProfs ! Number actually read

    ! Parameters

    CHARACTER (LEN=*), PARAMETER :: SZ_ERR = 'Failed to get size of &
         &dimension '
    CHARACTER (LEN=*), PARAMETER :: MLSMSG_INPUT = 'Error in input argument '
    CHARACTER (LEN=*), PARAMETER :: MLSMSG_L2GPRead = 'Unable to read L2GP &
                                                     &field:'
    ! Local Variables
    CHARACTER (LEN=80) :: list
    CHARACTER (LEN=480) :: msr

    INTEGER :: alloc_err, first, freq, lev, nDims, size, swid, status
    INTEGER :: start(3), stride(3), edge(3), dims(3)
    INTEGER :: nFreqs, nLevels, nTimes, nFreqsOr1, nLevelsOr1, myNumProfs

    LOGICAL :: firstCheck, lastCheck

    REAL, ALLOCATABLE :: realFreq(:), realSurf(:), realProf(:), real3(:,:,:)

    ! Attach to the swath for reading

    l2gp%Name = swathname

    swid = HE5_SWattach(L2FileHandle, l2gp%Name)
    IF (swid == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &attach to swath interface for reading.')

    ! Get dimension information

    lev = 0
    freq = 0

    nDims = HE5_SWinqdims(swid, list, dims)
    IF (nDims == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &get dimension information.')
    IF ( INDEX(list,'nLevels') /= 0 ) lev = 1
    IF ( INDEX(list,'Freq') /= 0 ) freq = 1

    size = HE5_SWdiminfo(swid, DIM_NAME1)
    IF (size == -1) THEN
       msr = SZ_ERR // DIM_NAME1
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%nTimes = size
    nTimes=size
    IF (lev == 0) THEN
       nLevels = 0
    ELSE
       size = HE5_SWdiminfo(swid, DIM_NAME2)
       IF (size == -1) THEN
          msr = SZ_ERR // DIM_NAME2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       nLevels = size

    ENDIF

    IF (freq == 1) THEN
       size = HE5_SWdiminfo(swid, DIM_NAME3)
       IF (size == -1) THEN
          msr = SZ_ERR // DIM_NAME3
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       nFreqs = size
    ELSE
       nFreqs = 0
    ENDIF

    ! Check optional input arguments

    firstCheck = PRESENT(firstProf)
    lastCheck = PRESENT(lastProf)

    IF (firstCheck) THEN

       IF ( (firstProf >= l2gp%nTimes) .OR. (firstProf < 0) ) THEN
          msr = MLSMSG_INPUT // 'firstProf'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ELSE
          first = firstProf
       ENDIF

    ELSE

       first = 0

    ENDIF

    IF (lastCheck) THEN

       IF (lastProf < first) THEN
          msr = MLSMSG_INPUT // 'lastProf'
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

       IF (lastProf >= nTimes) THEN
          myNumProfs = nTimes - first
       ELSE
          myNumProfs = lastProf - first + 1
       ENDIF

    ELSE

       myNumProfs = nTimes - first

    ENDIF

    ! Allocate result

    CALL SetupNewL2GPRecord (l2gp, nFreqs=nFreqs, nLevels=nLevels, nTimes=mynumProfs)

    ! Allocate temporary arrays

    nFreqsOr1=MAX(nFreqs,1)
    nLevelsOr1=MAX(nLevels, 1)
    ALLOCATE(realProf(myNumProfs), realSurf(l2gp%nLevels), &
         realFreq(l2gp%nFreqs), &
         real3(nFreqsOr1,nLevelsOr1,myNumProfs), STAT=alloc_err)

    ! Read the horizontal geolocation fields

    start(1) = 0
    start(2) = 0
    start(3) = first
    stride = 1
    edge(1) = nFreqsOr1
    edge(2) = nLevelsOr1
    edge(3) = myNumProfs
    print*,char(15)
!    print*,"Reading lats: start=",start
!    print*,"Stride=",stride
!    print*,"Edge=",edge
    status = HE5_SWrdfld(swid, GEO_FIELD1, start(3:3), stride(3:3), &
         edge(3:3), realProf)

    !print*,"Lats:",realProf( (/1,2,numProfs-1,numProfs /) )

    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD1
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%latitude = DBLE(realProf)

    status = HE5_SWrdfld(swid, GEO_FIELD2, start(3:3), stride(3:3), edge(3:3), &
         realProf)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD2
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%longitude = DBLE(realProf)

    status = HE5_SWrdfld(swid, GEO_FIELD3, start(3:3), stride(3:3), edge(3:3), &
         l2gp%time)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD3
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    status = HE5_SWrdfld(swid, GEO_FIELD4, start(3:3), stride(3:3), edge(3:3), &
         realProf)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD4
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    l2gp%solarTime = DBLE(realProf)
    status = HE5_SWrdfld(swid, GEO_FIELD5, start(3:3), stride(3:3), edge(3:3), &
         realProf)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD5
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%solarZenith = DBLE(realProf)

    status = HE5_SWrdfld(swid, GEO_FIELD6, start(3:3), stride(3:3), edge(3:3), &
         realProf)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD6
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%losAngle = DBLE(realProf)

    status = HE5_SWrdfld(swid, GEO_FIELD7, start(3:3), stride(3:3), edge(3:3), &
         realProf)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD7
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%geodAngle = DBLE(realProf)

    status = HE5_SWrdfld(swid, GEO_FIELD8, start(3:3), stride(3:3), edge(3:3), &
         l2gp%chunkNumber)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD8
       CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
    ENDIF

    ! Read the pressures vertical geolocation field, if it exists

    IF (lev /= 0) THEN

       status = HE5_SWrdfld(swid, GEO_FIELD9, start(2:2), stride(2:2), edge(2:2), &
            realSurf)
       IF (status == -1) THEN
          msr = MLSMSG_L2GPRead // GEO_FIELD9
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF

       l2gp%pressures = DBLE(realSurf)

    ENDIF

    ! Read the frequency geolocation field, if it exists

    IF (freq == 1) THEN

       edge(1) = l2gp%nFreqs

       status = HE5_SWrdfld(swid, GEO_FIELD10, start(1:1), stride(1:1), edge(1:1), &
            realFreq)
       IF (status == -1) THEN
          msr = MLSMSG_L2GPRead // GEO_FIELD10
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       l2gp%frequency = DBLE(realFreq)

    ENDIF

    ! Read the data fields that may have 1-3 dimensions

    IF ( freq == 1) THEN

       status = HE5_SWrdfld(swid, DATA_FIELD1, start, stride, edge, real3)
       IF (status == -1) THEN
          msr = MLSMSG_L2GPRead // DATA_FIELD1
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       l2gp%l2gpValue = DBLE(real3)

       status = HE5_SWrdfld(swid, DATA_FIELD2, start, stride, edge, real3)
       IF (status == -1) THEN
          msr = MLSMSG_L2GPRead // DATA_FIELD2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       l2gp%l2gpPrecision = DBLE(real3)

    ELSE IF ( lev == 1) THEN

       status = HE5_SWrdfld( swid, DATA_FIELD1, start(2:3), stride(2:3), &
            edge(2:3), real3 )
!            edge(2:3), real3(1,:,:) )
       IF (status == -1) THEN
          msr = MLSMSG_L2GPRead // DATA_FIELD1
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       l2gp%l2gpValue = DBLE(real3)

       status = HE5_SWrdfld( swid, DATA_FIELD2, start(2:3), stride(2:3), &
            edge(2:3), real3(1,:,:) )
       IF (status == -1) THEN
          msr = MLSMSG_L2GPRead // DATA_FIELD2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       l2gp%l2gpPrecision = DBLE(real3)

    ELSE

       status = HE5_SWrdfld( swid, DATA_FIELD1, start(3:3), stride(3:3), edge(3:3), &
            real3(1,1,:) )
       IF (status == -1) THEN
          msr = MLSMSG_L2GPRead // DATA_FIELD1
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       l2gp%l2gpValue = DBLE(real3)

       status = HE5_SWrdfld( swid, DATA_FIELD2, start(3:3), stride(3:3), edge(3:3), &
            real3(1,1,:) )
       IF (status == -1) THEN
          msr = MLSMSG_L2GPRead // DATA_FIELD2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       l2gp%l2gpPrecision = DBLE(real3)

    ENDIF
    
    ! Read the data fields that are 1-dimensional

    !         status = HE5_SWrdfld(swid, DATA_FIELD3,start(3:3),&
    !    stride(3:3),edge(3:3), l2gp%status)
    print*,"Warning: reading of status field disabled"
    status=0
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // DATA_FIELD3
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
         
    status = HE5_SWrdfld(swid, DATA_FIELD4, start(3:3), stride(3:3),&
         edge(3:3),realProf)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // DATA_FIELD4
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%quality = DBLE(realProf)

    ! Deallocate local variables

    DEALLOCATE(realFreq, realSurf, realProf, real3, STAT=alloc_err)
    IF ( alloc_err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
         'Failed deallocation of local real variables.')

    !  After reading, detach from HE5_SWath interface

    status = HE5_SWdetach(swid)
    IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &detach from swath interface after reading.')

    ! Set numProfs if wanted
    IF (PRESENT(numProfs)) numProfs=myNumProfs

    !-----------------------------
  END SUBROUTINE ReadL2GPData
  !-----------------------------

  ! --------------------------------------  OutputL2GP_createFile  -----
  SUBROUTINE OutputL2GP_createFile (l2gp, L2FileHandle, swathName,nLevels)

    ! Brief description of subroutine
    ! This subroutine sets up the structural definitions in an empty L2GP file.

    ! Arguments

    INTEGER, INTENT(in) :: L2FileHandle ! From swopen
    TYPE( L2GPData_T ), INTENT(inout) :: l2gp
    CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: swathName ! Defaults to l2gp%swathName
    INTEGER,optional::nLevels
    ! Parameters

    CHARACTER (len=*), PARAMETER :: DIM_ERR = 'Failed to define dimension '
    CHARACTER (len=*), PARAMETER :: GEO_ERR = &
         & 'Failed to define geolocation field '
    CHARACTER (len=*), PARAMETER :: DAT_ERR = 'Failed to define data field '

    ! Variables

    CHARACTER (len=480) :: MSR
    CHARACTER (len=132) :: NAME   ! From l2gp%name

    ! THESE ARE HDF5 CHUNKS, _NOT_ MLS ALONG-TRACK PROCESSING CHUNKS 
    INTEGER,DIMENSION(7)::CHUNK_DIMS
    INTEGER::CHUNK_RANK
    INTEGER::CHUNKTIMES,CHUNKFREQS,CHUNKLEVELS

    INTEGER :: SWID, STATUS
    character(len=1)::poop
    IF (PRESENT(swathName)) THEN
       name=swathName
    ELSE
       name=l2gp%name
    ENDIF
    chunktimes=1
    chunkfreqs=1 ! better as nFreqs, but I have yet to see a case with nfreqs>1
    if(present(nLevels))then
       chunklevels=nLevels
    else
       chunklevels=5
    endif
    
    ! Create the swath within the file
    print*,"Creating swath called ",name
    swid = HE5_SWcreate(L2FileHandle, TRIM(name))
    print*,"Swath ",name,"has SW id :",swid
    IF ( swid == -1 ) THEN
       msr = 'Failed to create swath ' // TRIM(name)
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    ! Define dimensions

    ! Defining special "unlimited dimension called UNLIM
    print*,"Defined Unlim with size", H5S_UNLIMITED
    status = HE5_SWdefdim(swid, UNLIM, H5S_UNLIMITED)

    print*,"Defining dimension ", DIM_NAME1," with size",l2gp%nTimes
    status = HE5_SWdefdim(swid, DIM_NAME1, l2gp%nTimes)


    IF ( status == -1 ) THEN
       msr = DIM_ERR // DIM_NAME1
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    IF ( l2gp%nLevels > 0 ) THEN
       
       print*,"Defining dimension ", DIM_NAME2," with size",l2gp%nLevels
       status = HE5_SWdefdim(swid, DIM_NAME2, l2gp%nLevels)
       IF ( status == -1 ) THEN
          msr = DIM_ERR // DIM_NAME2
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    IF ( l2gp%nFreqs > 0 ) THEN
       print*,"Defining dimension ", DIM_NAME3," with size",l2gp%nFreqs
       status = HE5_SWdefdim(swid, DIM_NAME3, l2gp%nFreqs)
       IF ( status == -1 ) THEN
          msr = DIM_ERR // DIM_NAME3
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    ! Define horizontal geolocation fields using above dimensions

    print*,"Defining geolocation field ",GEO_FIELD1," of dim. ", DIM_NAME1
    print*,"... and of type ",H5T_NATIVE_FLOAT
    chunk_rank=1
    chunk_dims=1
    chunk_dims(1)=CHUNKTIMES
    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    print*,"Set chunking -- status=",status
    status = HE5_SWdefgfld(swid, GEO_FIELD1, DIM_NAME1,MAX_DIML1,&
         H5T_NATIVE_FLOAT , 0)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD1
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF
    print*,"Defined geolocation field ",GEO_FIELD1,"of dim.", DIM_NAME1
    print*,"... and of type ",H5T_NATIVE_FLOAT

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)

    status = HE5_SWdefgfld(swid, GEO_FIELD2, DIM_NAME1, MAX_DIML1,&
         H5T_NATIVE_FLOAT, HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD2
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    status = HE5_SWdefgfld(swid, GEO_FIELD3, DIM_NAME1, MAX_DIML1, &
    H5T_NATIVE_DOUBLE, HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD3
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    status = HE5_SWdefgfld(swid, GEO_FIELD4, DIM_NAME1,MAX_DIML1,&
         H5T_NATIVE_FLOAT, HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD4
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    status = HE5_SWdefgfld(swid, GEO_FIELD5, DIM_NAME1, MAX_DIML1, &
         H5T_NATIVE_FLOAT, HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD5
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    status = HE5_SWdefgfld(swid, GEO_FIELD6, DIM_NAME1,MAX_DIML1,&
         H5T_NATIVE_FLOAT, HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD6
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    status = HE5_SWdefgfld(swid, GEO_FIELD7, DIM_NAME1, MAX_DIML1,&
    H5T_NATIVE_FLOAT,   HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD7
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    status = HE5_SWdefgfld(swid, GEO_FIELD8, DIM_NAME1, MAX_DIML1,&
         H5T_NATIVE_FLOAT, HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD8
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    IF ( l2gp%nLevels > 0 ) THEN

       status = HE5_SWdefgfld(swid, GEO_FIELD9, DIM_NAME2,MAX_DIML,&
            H5T_NATIVE_FLOAT, HDFE_NOMERGE)
       IF ( status == -1 ) THEN
          msr = GEO_ERR // GEO_FIELD9
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    IF ( l2gp%nFreqs > 0 ) THEN

       status = HE5_SWdefgfld(swid, GEO_FIELD10, DIM_NAME3,MAX_DIML,&
            H5T_NATIVE_FLOAT, HDFE_NOMERGE)
       IF ( status == -1 ) THEN
          msr = GEO_ERR // GEO_FIELD10
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    ! Define data fields using above dimensions

    IF ( (l2gp%nFreqs > 0) .AND. (l2gp%nLevels > 0) ) THEN
       chunk_rank=3
       chunk_dims(1:3)=(/ CHUNKFREQS,CHUNKLEVELS,CHUNKTIMES /)
       status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)

       status = HE5_SWdefdfld(swid, DATA_FIELD1, DIM_NAME123, MAX_DIML123,&
       H5T_NATIVE_FLOAT,HDFE_NOMERGE)

       IF ( status == -1 ) THEN
          msr = DAT_ERR // DATA_FIELD1 // ' for 3D quantity.'
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF


       status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
       status = HE5_SWdefdfld(swid, DATA_FIELD2, DIM_NAME123, MAX_DIML123,&
            H5T_NATIVE_FLOAT, HDFE_NOMERGE)

       IF ( status == -1 ) THEN
          msr = DAT_ERR // DATA_FIELD2 // ' for 3D quantity.'
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF


    ELSE IF ( l2gp%nLevels > 0 ) THEN
       chunk_rank=2
       chunk_dims(1:7)=(/ CHUNKLEVELS,CHUNKTIMES,37,38,39,47,49/)
       status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
       print*,"Set chunking with status=",status
       print*,"chunking=",chunk_dims
       print*,"About to define 2-D extendible field"
       !read(*,'(a1)')poop
       !print*,poop

       print*,"Calling SWdefdfld with args ",swid, DATA_FIELD1, &
            DIM_NAME12, MAX_DIML12, H5T_NATIVE_FLOAT, HDFE_NOMERGE
       status = HE5_SWdefdfld(swid, DATA_FIELD1, DIM_NAME12, MAX_DIML12, &
            H5T_NATIVE_FLOAT, HDFE_NOMERGE)
       print*,"Defined 2-D extendible field"
       !read(*,'(a1)')poop
       !print*,poop
    
       IF ( status == -1 ) THEN
          msr = DAT_ERR // DATA_FIELD1 //  ' for 2D quantity.'
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF

           status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
           status = HE5_SWdefdfld(swid, DATA_FIELD2, DIM_NAME12, MAX_DIML12,&
            H5T_NATIVE_FLOAT,HDFE_NOMERGE)

       IF ( status == -1 ) THEN
          msr = DAT_ERR // DATA_FIELD2 //  ' for 2D quantity.'
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF

    ELSE
       chunk_rank=1
       chunk_dims(1)=CHUNKTIMES
       status=HE5_SWdefchunk(swid,chunk_rank,chunk_dims)
    
       status = HE5_SWdefdfld(swid, DATA_FIELD1, DIM_NAME1,MAX_DIML1,&
            H5T_NATIVE_FLOAT,HDFE_NOMERGE)

       IF ( status == -1 ) THEN
          msr = DAT_ERR // DATA_FIELD1 // ' for 1D quantity.'
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF

       status = HE5_SWdefdfld(swid, DATA_FIELD2, DIM_NAME1, MAX_DIML1,&
       H5T_NATIVE_FLOAT, HDFE_NOMERGE)

       IF ( status == -1 ) THEN
          msr = DAT_ERR // DATA_FIELD2 // ' for 1D quantity.'
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF

    END IF

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
    IF ( status == -1 ) THEN
       msr = DAT_ERR // DATA_FIELD4
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    ! Detach from the HE5_SWath interface.This stores the swath info within the
    ! file and must be done before writing or reading data to or from the
    ! swath. (May be un-necessary for HDF5 -- test program works OK without.)
    ! 
    status = HE5_SWdetach(swid)
    IF ( status == -1 ) THEN
       CALL MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Failed to detach from swath interface after definition.' )
    END IF

    !--------------------------------------
  END SUBROUTINE OutputL2GP_createFile
  !--------------------------------------

  !-----------------------------------------  OutputL2GP_writeGeo  -----
  SUBROUTINE OutputL2GP_writeGeo (l2gp, l2FileHandle, swathName,offset)

    ! Brief description of subroutine
    ! This subroutine writes the geolocation fields to an L2GP output file.

    ! Arguments

    TYPE( L2GPData_T ), INTENT(inout) :: l2gp
    INTEGER, INTENT(in) :: l2FileHandle ! From swopen
    CHARACTER (len=*), INTENT(IN), OPTIONAL :: swathName ! Defaults->l2gp%name
    INTEGER,INTENT(IN),OPTIONAL::offset
    ! Parameters

    CHARACTER (len=*), PARAMETER :: WR_ERR = &
         & 'Failed to write geolocation field '
    
    ! Variables

    CHARACTER (len=480) :: msr
    CHARACTER (len=132) :: name ! Either swathName or l2gp%name
    
    INTEGER :: status, swid,myOffset
    INTEGER :: start(2), stride(2), edge(2)

    IF (PRESENT(offset)) THEN
       myOffset=offset
    ELSE
       myOffset=0
    ENDIF


    IF (PRESENT(swathName)) THEN
       name=swathName
    ELSE
       name=l2gp%name
    ENDIF

    swid = HE5_SWattach (l2FileHandle, name)

    ! Write data to the fields

    stride = 1
    start = myOffset
    edge(1) = l2gp%nTimes
    print*,"writeGeo Attached swath ",name," with SW ID=",swid
    print*,"About to write latitude with offset=",myoffset
    status = HE5_SWwrfld(swid, GEO_FIELD1, start, stride, edge, &
         REAL(l2gp%latitude))
    print*,"wrote latitude, maybe"
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD1
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = HE5_SWwrfld(swid, GEO_FIELD2, start, stride, edge, &
         REAL(l2gp%longitude))
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD2
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = HE5_SWwrfld(swid, GEO_FIELD3, start, stride, edge, &
         l2gp%time)
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD3
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF


    !print*,"writing ", REAL(l2gp%solarZenith)," as SZA"
    status = HE5_SWwrfld(swid, GEO_FIELD5, start, stride, edge, &
         REAL(l2gp%solarZenith))
    !print*,"just wrote ", REAL(l2gp%solarZenith)," as SZA"
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD5
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = HE5_SWwrfld(swid, GEO_FIELD4, start, stride, edge, &
        REAL(l2gp%solarTime))
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD4
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF


    status = HE5_SWwrfld(swid, GEO_FIELD6, start, stride, edge, &
         REAL(l2gp%losAngle))
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD6
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = HE5_SWwrfld(swid, GEO_FIELD7, start, stride, edge, &
         REAL(l2gp%geodAngle))
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD7
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = HE5_SWwrfld(swid, GEO_FIELD8, start, stride, edge, &
         l2gp%chunkNumber)
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD8
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    IF ( l2gp%nLevels > 0 ) THEN
       edge(1) = l2gp%nLevels
       start(1)=0 ! needed because offset may have made this /=0
       status = HE5_SWwrfld(swid, GEO_FIELD9, start, stride, edge, &
            REAL(l2gp%pressures))
       IF ( status == -1 ) THEN
          msr = WR_ERR // GEO_FIELD9
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    IF ( l2gp%nFreqs > 0 ) THEN
       edge(1) = l2gp%nFreqs
       start(1)=0 ! needed because offset may have made this /=0
       l2gp%frequency = 0
       status = HE5_SWwrfld(swid, GEO_FIELD10, start, stride, edge, &
            REAL(l2gp%frequency))
       IF ( status == -1 ) THEN
          msr = WR_ERR // GEO_FIELD10
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    ! Detach from the swath interface.  

    status = HE5_SWdetach(swid)
    print*,"Detatched from swath -- error=",status
    IF ( status == -1 ) THEN
       CALL MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach from swath interface' )
    END IF

    !------------------------------------
  END SUBROUTINE OutputL2GP_writeGeo
  !------------------------------------

  !----------------------------------------  OutputL2GP_writeData  -----
  SUBROUTINE OutputL2GP_writeData(l2gp, l2FileHandle, swathName,offset)

    ! Brief description of subroutine
    ! This subroutine writes the data fields to an L2GP output file.
    ! For now, you have to write all of l2gp, but you can choose to write
    ! it at some offset into the file
    ! Arguments

    TYPE( L2GPData_T ), INTENT(inout) :: l2gp
    INTEGER, INTENT(in) :: l2FileHandle ! From swopen
    CHARACTER (len=*), INTENT(IN), OPTIONAL :: swathName ! Defaults->l2gp%name
    INTEGER,INTENT(IN),OPTIONAL::offset
    ! Parameters

    CHARACTER (len=*), PARAMETER :: WR_ERR = 'Failed to write data field '

    ! Variables

    CHARACTER (len=480) :: msr
    CHARACTER (len=132) :: name     ! Either swathName or l2gp%name

    INTEGER :: status,myOffset
    INTEGER :: start(3), stride(3), edge(3)
    INTEGER :: swid
    
    IF (PRESENT(offset)) THEN
       myOffset=offset
    ELSE
       myOffset=0
    ENDIF

    IF (PRESENT(swathName)) THEN
       name=swathName
    ELSE
       name=l2gp%name
    ENDIF

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
    IF ( l2gp%nFreqs > 0 ) THEN
       !print*,"Writing 3D field"
       ! Value and Precision are 3-D fields
       status = HE5_SWwrfld(swid, DATA_FIELD1, start, stride, edge, &
            & RESHAPE(l2gp%l2gpValue, (/SIZE(l2gp%l2gpValue)/)) )
       IF ( status == -1 ) THEN
          msr = WR_ERR // DATA_FIELD1
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
       status = HE5_SWwrfld(swid, DATA_FIELD2, start, stride, edge, &
            & RESHAPE(REAL(l2gp%l2gpPrecision), (/SIZE(l2gp%l2gpPrecision)/)) )
       IF ( status == -1 ) THEN
          msr = WR_ERR // DATA_FIELD2
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF

    ELSE IF ( l2gp%nLevels > 0 ) THEN
       !Print*,"Writing 2-d field"
       ! Value and Precision are 2-D fields
      print*,"About to write data with offset=",myOffset
      
       status = HE5_SWwrfld( swid, DATA_FIELD1, start(2:3), stride(2:3), &
            edge(2:3), REAL(l2gp%l2gpValue(1,:,:) ))
       !print*,"Status of write was ",status
       IF ( status == -1 ) THEN
          msr = WR_ERR // DATA_FIELD1
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
       status = HE5_SWwrfld( swid, DATA_FIELD2, start(2:3), stride(2:3), &
            edge(2:3), REAL(l2gp%l2gpPrecision(1,:,:) ))
       IF ( status == -1 ) THEN
          msr = WR_ERR // DATA_FIELD2
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    ELSE

       ! Value and Precision are 1-D fields
       !Print*,"Writing 1-D field"
       status = HE5_SWwrfld( swid, DATA_FIELD1, start(3:3), stride(3:3), edge(3:3), &
            REAL(l2gp%l2gpValue(1,1,:) ))
       IF ( status == -1 ) THEN
          msr = WR_ERR // DATA_FIELD1
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
       status = HE5_SWwrfld( swid, DATA_FIELD2, start(3:3), stride(3:3), edge(3:3), &
            REAL(l2gp%l2gpPrecision(1,1,:) ))
       IF ( status == -1 ) THEN
          msr = WR_ERR // DATA_FIELD2
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    ! 1-D status & quality fields

    !HDF-EOS5 won't write a dataset of chars from FORTRAN
    !    status = HE5_SWwrfld(swid, DATA_FIELD3, start(3:3), stride(3:3),&
    !        edge(3:3), l2gp%status) ! 
    status=0
    print*,"Warning. Writing of status field disabled"
    IF ( status == -1 ) THEN
       msr = WR_ERR // DATA_FIELD3
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF
    l2gp%quality = 0
    status = HE5_SWwrfld(swid, DATA_FIELD4, start(3:3), stride(3:3), edge(3:3), &
         REAL(l2gp%quality))
    IF ( status == -1 ) THEN
       msr = WR_ERR // DATA_FIELD4
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    !     Detach from the swath interface.

    status = HE5_SWdetach(swid)
    print*,"Detatched from swath -- error=",status
    IF ( status == -1 ) THEN
       CALL MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach  from swath interface' )
    END IF


    !-------------------------------------
  END SUBROUTINE OutputL2GP_writeData
  !-------------------------------------

  ! --------------------------------------------------------------------------

  ! This subroutine is an amalgamation of the last three
  ! Should be renamed CreateAndWriteL2GPData
  SUBROUTINE WriteL2GPData(l2gp,l2FileHandle,swathName)

    ! Arguments

    INTEGER, INTENT(IN) :: l2FileHandle ! From swopen
    TYPE (L2GPData_T), INTENT(INOUT) :: l2gp
    CHARACTER (LEN=*), OPTIONAL, INTENT(IN) ::swathName!default->l2gp%swathName
    ! Exectuable code

    CALL OutputL2GP_createFile (l2gp, l2FileHandle, swathName)
    CALL OutputL2GP_writeGeo (l2gp, l2FileHandle, swathName)
    CALL OutputL2GP_writeData (l2gp, l2FileHandle, swathName)

  END SUBROUTINE WriteL2GPData
  !-------------------------------------------------------------


  SUBROUTINE AppendL2GPData(l2gp,l2FileHandle,swathName,offset)
    ! sticks l2gp into the swath swathName in the file pointed at by
    ! l2FileHandle,starting at the profile number "offset" (First profile
    ! in the file has offset==0). If this runs off the end ofthe swath, 
    ! it is lengthend automagically. 
    ! Arguments

    INTEGER, INTENT(IN) :: l2FileHandle ! From swopen
    TYPE (L2GPData_T), INTENT(INOUT) :: l2gp
    CHARACTER (LEN=*), OPTIONAL, INTENT(IN) ::swathName!default->l2gp%swathName
    INTEGER,INTENT(IN),OPTIONAL::offset
    !----Local variable
    INTEGER::myOffset
    ! ----Executable code---
    IF (PRESENT(offset)) THEN
       myOffset=offset
    ELSE
       myOffset=0
    ENDIF
    CALL OutputL2GP_writeGeo (l2gp, l2FileHandle, swathName,myOffset)
    CALL OutputL2GP_writeData (l2gp, l2FileHandle, swathName,myOffset)

  END SUBROUTINE AppendL2GPData



  ! ------------------------------------------ DUMP_L2GP ------------

  subroutine Dump_L2GP ( L2gp, Name )

    ! Dummy arguments
    type (l2gpData_T), intent(in) :: L2GP(:)
    character(len=*), optional :: Name

    ! Local variables
    integer :: i

    if ( present(name) ) call output ( name, advance='yes' )
    do i = 1, size(l2gp)
      call output ( 'L2GP Data: ')
      call display_string ( l2gp(i)%nameIndex, advance='yes' )
      call output ( 'nTimes: ')
      call output ( l2gp(i)%nTimes, 5)
      call output ( '  nLevels: ')
      call output ( l2gp(i)%nLevels, 3)
      call output ( '  nFreqs: ')
      call output ( l2gp(i)%nFreqs, 3, advance='yes')
      
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
END MODULE L2GPData
!=============================================================================

!
! $Log$
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
