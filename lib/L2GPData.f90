! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L2GPData                 ! Creation, manipulation and I/O for L2GP Data
  !=============================================================================

  USE Allocate_Deallocate, ONLY: Allocate_test, Deallocate_test
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
       & MLSMSG_Error, MLSMSG_Warning
  USE MLSCommon, ONLY: R8

  USE Hdf, ONLY: DFNT_CHAR8, DFNT_FLOAT32, DFNT_INT32, DFNT_FLOAT64
  USE HDFEOS!, ONLY: SWATTACH, SWCREATE, SWDEFDFLD, SWDEFDIM, SWDEFGFLD, &
     !& SWDETACH
  USE SWAPI, ONLY: SWWRFLD, SWRDFLD
  
  IMPLICIT NONE
!  INTEGER :: SWRDFLD
!  EXTERNAL SWRDFLD !Should USE SWAPI
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER(len=256), PRIVATE :: Id = &
       & "$Id$"
  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName = &
       & "$RCSfile$"
  !---------------------------------------------------------------------------


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
   CHARACTER (len=*), PARAMETER :: GEO_FIELD10= 'Frequency'

   CHARACTER (len=*), PARAMETER :: DIM_NAME1 = 'nTimes'
   CHARACTER (len=*), PARAMETER :: DIM_NAME2 = 'nLevels'
   CHARACTER (len=*), PARAMETER :: DIM_NAME3 = 'nFreqs'
   CHARACTER (len=*), PARAMETER :: DIM_NAME12 = 'nLevels,nTimes' ! In Fortran order?!!
   CHARACTER (len=*), PARAMETER :: DIM_NAME123 = 'nFreqs,nLevels,nTimes' ! as above

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
     REAL (r8), POINTER, DIMENSION(:) :: pressures ! Vertical coords (nLevels)

     ! Now the horizontal geolocation information. Dimensioned (nTimes)
     REAL (r8), POINTER, DIMENSION(:) :: latitude, longitude, solarTime, &
          & solarZenith, losAngle, geodAngle
     REAL (r8), POINTER, DIMENSION(:) :: time
     INTEGER, POINTER, DIMENSION(:) :: chunkNumber !

     ! Now we store the `frequency' geolocation field

     REAL (r8), POINTER, DIMENSION(:) :: frequency
     !        dimensioned (nFreqs)

     ! Finally we store the data fields

     REAL (r8), POINTER, DIMENSION(:,:,:) :: l2gpValue
     REAL (r8), POINTER, DIMENSION(:,:,:) :: l2gpPrecision
     ! dimensioned (nFreqs, nLevels, nTimes)

     CHARACTER (len=1), POINTER, DIMENSION(:) :: status
     REAL (r8), POINTER, DIMENSION(:) :: quality
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

    CALL SetupNewL2GPRecord( l2gp, nFreqs=l2gp%nFreqs, nLevels=l2gp%nLevels, nTimes=newNTimes)

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
    INTEGER, INTENT(OUT), OPTIONAL :: numProfs ! Number actually read

    ! Local Parameters
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

    swid = swattach(L2FileHandle, TRIM(l2gp%Name))
    IF (swid == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &attach to swath interface for reading.')

    ! Get dimension information

    lev = 0
    freq = 0

    nDims = swinqdims(swid, list, dims)
    IF (nDims == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &get dimension information.')
    IF ( INDEX(list,'nLevels') /= 0 ) lev = 1
    IF ( INDEX(list,'Freq') /= 0 ) freq = 1

    size = swdiminfo(swid, DIM_NAME1)
    IF (size == -1) THEN
       msr = SZ_ERR // DIM_NAME1
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%nTimes = size
    nTimes=size
    IF (lev == 0) THEN
       nLevels = 0
    ELSE
       size = swdiminfo(swid, DIM_NAME2)
       IF (size == -1) THEN
          msr = SZ_ERR // DIM_NAME2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       nLevels = size

    ENDIF

    IF (freq == 1) THEN
       size = swdiminfo(swid, DIM_NAME3)
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

    CALL SetupNewL2GPRecord (l2gp, nFreqs=nFreqs, nLevels=nLevels, nTimes=myNumProfs)

    ! Allocate temporary arrays

    nFreqsOr1=MAX(nFreqs,1)
    nLevelsOr1=MAX(nLevels, 1)
    ALLOCATE(realProf(myNumProfs), realSurf(l2gp%nLevels), &
         realFreq(l2gp%nFreqs), &
         real3(nFreqsOr1,nLevelsOr1,myNumProfs), STAT=alloc_err)
    IF (alloc_err /= 0) call MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
      & ' various things in ReadL2GPData')

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
          realProf)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD1
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%latitude = DBLE(realProf)

     status = swrdfld(swid, GEO_FIELD2, start(3:3), stride(3:3), edge(3:3), &
          realProf)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD2
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%longitude = DBLE(realProf)

     status = swrdfld(swid, GEO_FIELD3, start(3:3), stride(3:3), edge(3:3), &
          l2gp%time)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD3
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

     status = swrdfld(swid, GEO_FIELD4, start(3:3), stride(3:3), edge(3:3), &
          realProf)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD4
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%solarTime = DBLE(realProf)

     status = swrdfld(swid, GEO_FIELD5, start(3:3), stride(3:3), edge(3:3), &
          realProf)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD5
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%solarZenith = DBLE(realProf)

     status = swrdfld(swid, GEO_FIELD6, start(3:3), stride(3:3), edge(3:3), &
          realProf)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD6
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%losAngle = DBLE(realProf)

     status = swrdfld(swid, GEO_FIELD7, start(3:3), stride(3:3), edge(3:3), &
          realProf)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD7
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%geodAngle = DBLE(realProf)

     status = swrdfld(swid, GEO_FIELD8, start(3:3), stride(3:3), edge(3:3), &
          l2gp%chunkNumber)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // GEO_FIELD8
       CALL MLSMessage(MLSMSG_Warning, ModuleName, msr)
    ENDIF

    ! Read the pressures vertical geolocation field, if it exists

    IF (lev /= 0) THEN

       status = swrdfld(swid, GEO_FIELD9, start(2:2), stride(2:2), edge(2:2), &
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

       status = swrdfld(swid, GEO_FIELD10, start(1:1), stride(1:1), edge(1:1), &
         realFreq)
       IF (status == -1) THEN
          msr = MLSMSG_L2GPRead // GEO_FIELD10
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       l2gp%frequency = DBLE(realFreq)

    ENDIF

    ! Read the data fields that may have 1-3 dimensions

    IF ( freq == 1) THEN

       status = swrdfld(swid, DATA_FIELD1, start, stride, edge, real3)
       IF (status == -1) THEN
          msr = MLSMSG_L2GPRead // DATA_FIELD1
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       l2gp%l2gpValue = DBLE(real3)

       status = swrdfld(swid, DATA_FIELD2, start, stride, edge, real3)
       IF (status == -1) THEN
          msr = MLSMSG_L2GPRead // DATA_FIELD2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       l2gp%l2gpPrecision = DBLE(real3)

    ELSE IF ( lev == 1) THEN

      status = swrdfld( swid, DATA_FIELD1, start(2:3), stride(2:3), &
           edge(2:3), real3(1,:,:) )
      IF (status == -1) THEN
        msr = MLSMSG_L2GPRead // DATA_FIELD1
        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l2gp%l2gpValue = DBLE(real3)
      
      status = swrdfld( swid, DATA_FIELD2, start(2:3), stride(2:3), &
        edge(2:3), real3(1,:,:) )
      IF (status == -1) THEN
        msr = MLSMSG_L2GPRead // DATA_FIELD2
        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF
      l2gp%l2gpPrecision = DBLE(real3)
      
    ELSE

       status = swrdfld( swid, DATA_FIELD1, start(3:3), stride(3:3), edge(3:3), &
            real3(1,1,:) )
       IF (status == -1) THEN
          msr = MLSMSG_L2GPRead // DATA_FIELD1
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       l2gp%l2gpValue = DBLE(real3)

       status = swrdfld( swid, DATA_FIELD2, start(3:3), stride(3:3), edge(3:3), &
            real3(1,1,:) )
       IF (status == -1) THEN
          msr = MLSMSG_L2GPRead // DATA_FIELD2
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
       ENDIF
       l2gp%l2gpPrecision = DBLE(real3)

    ENDIF

    ! Read the data fields that are 1-dimensional

!    status = swrdfld(swid, DATA_FIELD3,start(3:3),stride(3:3),edge(3:3),&
!      l2gp%status)
! These lines commented out as they make NAG core dump on the deallocate statement.
! below.
    IF (status == -1) THEN
      msr = MLSMSG_L2GPRead // DATA_FIELD3
      CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF

    status = swrdfld(swid, DATA_FIELD4, start(3:3), stride(3:3), edge(3:3), &
         realProf)
    IF (status == -1) THEN
       msr = MLSMSG_L2GPRead // DATA_FIELD4
       CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
    ENDIF
    l2gp%quality = DBLE(realProf)


    ! Deallocate local variables

    DEALLOCATE(realFreq, realSurf, realProf, real3, STAT=alloc_err)
    IF ( alloc_err /= 0 ) CALL MLSMessage(MLSMSG_Error, ModuleName, &
         'Failed deallocation of local real variables.')

    !  After reading, detach from swath interface

    status = swdetach(swid)
    IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &detach from swath interface after reading.')

    ! Set numProfs if wanted
    IF (PRESENT(numProfs)) numProfs=myNumProfs

    !-----------------------------
  END SUBROUTINE ReadL2GPData
  !-----------------------------

  ! --------------------------------------  OutputL2GP_createFile  -----
  SUBROUTINE OutputL2GP_createFile (l2gp, L2FileHandle, swathName)

    ! Brief description of subroutine
    ! This subroutine sets up the structural definitions in an empty L2GP file.

    ! Arguments

    INTEGER, INTENT(in) :: L2FileHandle ! From swopen
    TYPE( L2GPData_T ), INTENT(inout) :: l2gp
    CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: swathName ! Defaults to l2gp%swathName

    ! Parameters

    CHARACTER (len=*), PARAMETER :: DIM_ERR = 'Failed to define dimension '
    CHARACTER (len=*), PARAMETER :: GEO_ERR = &
         & 'Failed to define geolocation field '
    CHARACTER (len=*), PARAMETER :: DAT_ERR = 'Failed to define data field '

    ! Variables

    CHARACTER (len=480) :: MSR
    CHARACTER (len=132) :: NAME   ! From l2gp%name

    INTEGER :: SWID, STATUS

    IF (PRESENT(swathName)) THEN
       name=swathName
    ELSE
       name=l2gp%name
    ENDIF

    ! Create the swath within the file

    swid = swcreate(L2FileHandle, TRIM(name))
    IF ( swid == -1 ) THEN
       msr = 'Failed to create swath ' // TRIM(name)
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    ! Define dimensions

    status = swdefdim(swid, DIM_NAME1, l2gp%nTimes)
    IF ( status == -1 ) THEN
       msr = DIM_ERR // DIM_NAME1
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    IF ( l2gp%nLevels > 0 ) THEN
       status = swdefdim(swid, DIM_NAME2, l2gp%nLevels)
       IF ( status == -1 ) THEN
          msr = DIM_ERR // DIM_NAME2
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    IF ( l2gp%nFreqs > 1 ) THEN
       status = swdefdim(swid, DIM_NAME3, l2gp%nFreqs)
       IF ( status == -1 ) THEN
          msr = DIM_ERR // DIM_NAME3
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    ! Define horizontal geolocation fields using above dimensions

    status = swdefgfld(swid, GEO_FIELD1, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD1
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swdefgfld(swid, GEO_FIELD2, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD2
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swdefgfld(swid, GEO_FIELD3, DIM_NAME1, DFNT_FLOAT64, &
         HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD3
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swdefgfld(swid, GEO_FIELD4, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD4
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swdefgfld(swid, GEO_FIELD5, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD5
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swdefgfld(swid, GEO_FIELD6, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD6
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swdefgfld(swid, GEO_FIELD7, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD7
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swdefgfld(swid, GEO_FIELD8, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = GEO_ERR // GEO_FIELD8
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    IF ( l2gp%nLevels > 0 ) THEN
       status = swdefgfld(swid, GEO_FIELD9, DIM_NAME2, DFNT_FLOAT32, &
            HDFE_NOMERGE)
       IF ( status == -1 ) THEN
          msr = GEO_ERR // GEO_FIELD9
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    IF ( l2gp%nFreqs > 0 ) THEN
       status = swdefgfld(swid, GEO_FIELD10, DIM_NAME3, DFNT_FLOAT32, &
            HDFE_NOMERGE)
       IF ( status == -1 ) THEN
          msr = GEO_ERR // GEO_FIELD10
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    ! Define data fields using above dimensions

    IF ( (l2gp%nFreqs > 0) .AND. (l2gp%nLevels > 0) ) THEN

       status = swdefdfld(swid, DATA_FIELD1, DIM_NAME123, DFNT_FLOAT32, &
            HDFE_NOMERGE)

       IF ( status == -1 ) THEN
          msr = DAT_ERR // DATA_FIELD1 // ' for 3D quantity.'
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF


       status = swdefdfld(swid, DATA_FIELD2, DIM_NAME123, DFNT_FLOAT32, &
            HDFE_NOMERGE)

       IF ( status == -1 ) THEN
          msr = DAT_ERR // DATA_FIELD2 // ' for 3D quantity.'
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF


    ELSE IF ( l2gp%nLevels > 0 ) THEN

       status = swdefdfld(swid, DATA_FIELD1, DIM_NAME12, DFNT_FLOAT32, &
            HDFE_NOMERGE)

       IF ( status == -1 ) THEN
          msr = DAT_ERR // DATA_FIELD1 //  ' for 2D quantity.'
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF

       status = swdefdfld(swid, DATA_FIELD2, DIM_NAME12, DFNT_FLOAT32, &
            HDFE_NOMERGE)

       IF ( status == -1 ) THEN
          msr = DAT_ERR // DATA_FIELD2 //  ' for 2D quantity.'
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF

    ELSE

       status = swdefdfld(swid, DATA_FIELD1, DIM_NAME1, DFNT_FLOAT32, &
            HDFE_NOMERGE)

       IF ( status == -1 ) THEN
          msr = DAT_ERR // DATA_FIELD1 // ' for 1D quantity.'
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF

       status = swdefdfld(swid, DATA_FIELD2, DIM_NAME1, DFNT_FLOAT32, &
            HDFE_NOMERGE)

       IF ( status == -1 ) THEN
          msr = DAT_ERR // DATA_FIELD2 // ' for 1D quantity.'
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF

    END IF

    status = swdefdfld(swid, DATA_FIELD3, DIM_NAME1, DFNT_CHAR8, &
         HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = DAT_ERR // DATA_FIELD3
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swdefdfld(swid, DATA_FIELD4, DIM_NAME1, DFNT_FLOAT32, &
         HDFE_NOMERGE)
    IF ( status == -1 ) THEN
       msr = DAT_ERR // DATA_FIELD4
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    ! Detach from the swath interface.  This stores the swath info within the
    ! file and must be done before writing or reading data to or from the
    ! swath.

    status = swdetach(swid)
    IF ( status == -1 ) THEN
       CALL MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Failed to detach from swath interface after definition.' )
    END IF

    !--------------------------------------
  END SUBROUTINE OutputL2GP_createFile
  !--------------------------------------

  !-----------------------------------------  OutputL2GP_writeGeo  -----
  SUBROUTINE OutputL2GP_writeGeo (l2gp, l2FileHandle, swathName)

    ! Brief description of subroutine
    ! This subroutine writes the geolocation fields to an L2GP output file.

    ! Arguments

    TYPE( L2GPData_T ), INTENT(inout) :: l2gp
    INTEGER, INTENT(in) :: l2FileHandle ! From swopen
    CHARACTER (len=*), INTENT(IN), OPTIONAL :: swathName ! Defaults to l2gp%name

    ! Parameters

    CHARACTER (len=*), PARAMETER :: WR_ERR = &
         & 'Failed to write geolocation field '

    ! Variables

    CHARACTER (len=480) :: msr
    CHARACTER (len=132) :: name ! Either swathName or l2gp%name

    INTEGER :: status, swid
    INTEGER :: start(2), stride(2), edge(2)

    IF (PRESENT(swathName)) THEN
       name=swathName
    ELSE
       name=l2gp%name
    ENDIF

    swid = swattach (l2FileHandle, name)

    ! Write data to the fields

    stride(1) = 1
    start(1) = 0
    edge(1) = l2gp%nTimes

    status = swwrfld(swid, GEO_FIELD1, start, stride, edge, &
         REAL(l2gp%latitude))
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD1
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swwrfld(swid, GEO_FIELD2, start, stride, edge, &
         REAL(l2gp%longitude))
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD2
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swwrfld(swid, GEO_FIELD3, start, stride, edge, &
         l2gp%time)
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD3
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swwrfld(swid, GEO_FIELD4, start, stride, edge, &
        REAL(l2gp%solarTime))
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD4
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swwrfld(swid, GEO_FIELD5, start, stride, edge, &
         REAL(l2gp%solarZenith))
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD5
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swwrfld(swid, GEO_FIELD6, start, stride, edge, &
         REAL(l2gp%losAngle))
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD6
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swwrfld(swid, GEO_FIELD7, start, stride, edge, &
         REAL(l2gp%geodAngle))
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD7
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    status = swwrfld(swid, GEO_FIELD8, start, stride, edge, &
         l2gp%chunkNumber)
    IF ( status == -1 ) THEN
       msr = WR_ERR // GEO_FIELD8
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    IF ( l2gp%nLevels > 0 ) THEN
       edge(1) = l2gp%nLevels
       status = swwrfld(swid, GEO_FIELD9, start, stride, edge, &
            REAL(l2gp%pressures))
       IF ( status == -1 ) THEN
          msr = WR_ERR // GEO_FIELD9
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    IF ( l2gp%nFreqs > 0 ) THEN
       edge(1) = l2gp%nFreqs
       l2gp%frequency = 0
       status = swwrfld(swid, GEO_FIELD10, start, stride, edge, &
            REAL(l2gp%frequency))
       IF ( status == -1 ) THEN
          msr = WR_ERR // GEO_FIELD10
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    ! Detach from the swath interface.  

    status = swdetach(swid)

    IF ( status == -1 ) THEN
       CALL MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach from swath interface' )
    END IF

    !------------------------------------
  END SUBROUTINE OutputL2GP_writeGeo
  !------------------------------------

  !----------------------------------------  OutputL2GP_writeData  -----
  SUBROUTINE OutputL2GP_writeData(l2gp, l2FileHandle, swathName)

    ! Brief description of subroutine
    ! This subroutine writes the data fields to an L2GP output file.

    ! Arguments

    TYPE( L2GPData_T ), INTENT(inout) :: l2gp
    INTEGER, INTENT(in) :: l2FileHandle ! From swopen
    CHARACTER (len=*), INTENT(IN), OPTIONAL :: swathName ! Defaults to l2gp%name

    ! Parameters

    CHARACTER (len=*), PARAMETER :: WR_ERR = 'Failed to write data field '

    ! Variables

    CHARACTER (len=480) :: msr
    CHARACTER (len=132) :: name     ! Either swathName or l2gp%name

    INTEGER :: status
    INTEGER :: start(3), stride(3), edge(3)
    INTEGER :: swid

    IF (PRESENT(swathName)) THEN
       name=swathName
    ELSE
       name=l2gp%name
    ENDIF
    ! Write data to the fields

    start = 0
    stride = 1
    edge(1) = l2gp%nFreqs
    edge(2) = l2gp%nLevels
    edge(3) = l2gp%nTimes
    swid = swattach (l2FileHandle, name)
    IF ( l2gp%nFreqs > 0 ) THEN
       ! Value and Precision are 3-D fields

       status = swwrfld(swid, DATA_FIELD1, start, stride, edge, &
            & RESHAPE(l2gp%l2gpValue, (/SIZE(l2gp%l2gpValue)/)) )
       IF ( status == -1 ) THEN
          msr = WR_ERR // DATA_FIELD1
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
       status = swwrfld(swid, DATA_FIELD2, start, stride, edge, &
            & RESHAPE(REAL(l2gp%l2gpPrecision), (/SIZE(l2gp%l2gpPrecision)/)) )
       IF ( status == -1 ) THEN
          msr = WR_ERR // DATA_FIELD2
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF

    ELSE IF ( l2gp%nLevels > 0 ) THEN
       ! Value and Precision are 2-D fields

       status = swwrfld( swid, DATA_FIELD1, start(2:3), stride(2:3), &
            edge(2:3), REAL(l2gp%l2gpValue(1,:,:)) )

       IF ( status == -1 ) THEN
          msr = WR_ERR // DATA_FIELD1
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
       status = swwrfld( swid, DATA_FIELD2, start(2:3), stride(2:3), &
            edge(2:3), REAL(l2gp%l2gpPrecision(1,:,:) ))
       IF ( status == -1 ) THEN
          msr = WR_ERR // DATA_FIELD2
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    ELSE

       ! Value and Precision are 1-D fields
       status = swwrfld( swid, DATA_FIELD1, start(3:3), stride(3:3), edge(3:3), &
            REAL(l2gp%l2gpValue(1,1,:) ))
       IF ( status == -1 ) THEN
          msr = WR_ERR // DATA_FIELD1
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
       status = swwrfld( swid, DATA_FIELD2, start(3:3), stride(3:3), edge(3:3), &
            REAL(l2gp%l2gpPrecision(1,1,:) ))
       IF ( status == -1 ) THEN
          msr = WR_ERR // DATA_FIELD2
          CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
       END IF
    END IF

    ! 1-D status & quality fields

    status = swwrfld(swid, DATA_FIELD3, start(3:3), stride(3:3), edge(3:3), &
         l2gp%status)  
    IF ( status == -1 ) THEN
       msr = WR_ERR // DATA_FIELD3
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF
    l2gp%quality = 0
    status = swwrfld(swid, DATA_FIELD4, start(3:3), stride(3:3), edge(3:3), &
         REAL(l2gp%quality))
    IF ( status == -1 ) THEN
       msr = WR_ERR // DATA_FIELD4
       CALL MLSMessage ( MLSMSG_Error, ModuleName, msr )
    END IF

    !     Detach from the swath interface.

    status = swdetach(swid)
    IF ( status == -1 ) THEN
       CALL MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach  from swath interface' )
    END IF


    !-------------------------------------
  END SUBROUTINE OutputL2GP_writeData
  !-------------------------------------

  ! --------------------------------------------------------------------------

  ! This subroutine is an amalgamation of the last three

  SUBROUTINE WriteL2GPData(l2gp,l2FileHandle,swathName)

    ! Arguments

    INTEGER, INTENT(IN) :: l2FileHandle ! From swopen
    TYPE (L2GPData_T), INTENT(INOUT) :: l2gp
    CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: swathName ! (defaults to l2gp%swathName)

    ! Exectuable code

    CALL OutputL2GP_createFile (l2gp, l2FileHandle, swathName)
    CALL OutputL2GP_writeGeo (l2gp, l2FileHandle, swathName)
    CALL OutputL2GP_writeData (l2gp, l2FileHandle, swathName)

  END SUBROUTINE WriteL2GPData

  !=============================================================================
END MODULE L2GPData
!=============================================================================

!
! $Log$
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
