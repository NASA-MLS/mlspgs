! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE SignalsFile              ! Read emls-signals.dat file
!============================================================================
  
  USE MLSStrings
  
  IMPLICIT NONE
  PUBLIC
  
  PRIVATE :: id
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = & 
     "$Id$"
  !-----------------------------------------------------------------------------
  
  
  ! This module reads and parses the file describing the valid signals for the
  ! MLS instrument, currently this is called emls-signals.dat, but the name may
  ! change of course.
  
  ! First we define the global parameters and data types, we'll prefix these
  ! all with SDB for Signals Database.

  INTEGER, PARAMETER :: SDB_NameLen=10

  ! This datatype describes a radiometer
  TYPE SDB_RadiometerInfo_T
     CHARACTER (LEN=SDB_NameLen) :: name   ! e.g. R1A:118
     CHARACTER (LEN=SDB_NameLen) :: prefix ! e.g. R1A
     CHARACTER (LEN=SDB_NameLen) :: suffix ! e.g. 118
     INTEGER :: number                     ! e.g. 1
     CHARACTER(LEN=1) :: modifier          ! e.g. A/B or H/V for R5 (emls)
     DOUBLE PRECISION :: lo                ! Local oscillator /MHz
  END TYPE SDB_RadiometerInfo_T

  ! This datatype describes a band
  TYPE SDB_BandInfo_T
     CHARACTER (LEN=SDB_NameLen) :: name     ! e.g. B1F:PT
     CHARACTER (LEN=SDB_NameLen) :: suffix   ! e.g. PT
     INTEGER :: number                       ! e.g. 1
     CHARACTER (LEN=1) :: spectrometerFamily ! e.g. F
     INTEGER :: spectrometerFamilyIndex      ! Index into array of next type
     INTEGER :: specificChannelIndex         ! Index into array for WF4 chans
     DOUBLE PRECISION :: centerFreqIF        ! Center i.f. frequency (MHz)
  END TYPE SDB_BandInfo_T

  ! This datatype describes a spectrometer family
  TYPE SDB_SpectrometerFamilyInfo_T
     CHARACTER (LEN=SDB_NameLen) :: name ! Name of family e.g. FB25
     INTEGER :: noSpectrometers ! Number of spectrometers in this family
     INTEGER :: noChannels      ! Number of channels in family
     INTEGER :: firstChannel    ! First channel number (e.g. 1)
     INTEGER :: lastChannel     ! Last channel number (e.g. 25)
     LOGICAL :: specific        ! If set have discrete freqs. e.g. wf4 series

     ! These two arrays are the position and width of the channels wrt. the if
     ! The arrays are actually dimensioned firstChannel:lastChannel.
     ! Units are MHz

     DOUBLE PRECISION, DIMENSION(:), POINTER :: position
     DOUBLE PRECISION, DIMENSION(:), POINTER :: width
  END TYPE SDB_SpectrometerFamilyInfo_T

  ! This datatype describes a valid radiometer/switch/band/spectrometer
  ! combination
  TYPE SDB_ValidSignal_T
     INTEGER :: radiometerIndex ! Index into array of radiometerInfos
     INTEGER :: bandIndex       ! Index into array of bandInfos
     INTEGER :: switchIndex     ! Index into array of switch names
     INTEGER :: spectrometerFamilyIndex ! Index into array of spec. fams.
     INTEGER :: spectrometerNumber ! Note count from *one* for this!
  END TYPE SDB_ValidSignal_T

  ! This datatype is an amalgam of the above and is the database that is filled
  TYPE MLS_SignalsDatabase_T
     INTEGER :: noRadiometers   ! Including redundant etc.
     INTEGER :: noBands         ! Accross the whole instrument
     INTEGER :: noSwitches      ! No. vald S0, S1 etc. fields
     INTEGER :: noSpectrometerFamilies
     INTEGER :: noValidSignals    ! Number of valid Signal combinations
     INTEGER :: noSpecificChannels ! Total no. channels in WF4 type specs.
     
     TYPE (SDB_RadiometerInfo_T), DIMENSION(:), POINTER :: radiometerInfo
                                ! Actually dimensioned (0:noRadiometers-1)
     TYPE (SDB_BandInfo_T), DIMENSION(:), POINTER :: bandIndo
                                ! Actually dimensioned (0:noBands-1)
     CHARACTER (LEN=SDB_NameLen), DIMENSION(:), POINTER :: switches
                                ! Actually dimensioned (0:noSwitches-1)
     TYPE (SDB_SpectrometerFamilyInfo_T), DIMENSION(:), POINTER :: &
          & spectrometerFamilies
             ! Actually dimensioned (0:noSpectrometerFamilies-1)
     TYPE (SDB_ValidSignal_T), DIMENSION(:), POINTER :: validSignal
                                ! Valid Signal combinations
                                ! Actually dimensioned (0:noValidSignals-1)

     ! These two are for `specific' channels, they describe their position and
     ! width in i.f. space in MHz.
     DOUBLE PRECISION, DIMENSION(:), POINTER :: specificPosition, specificWidth
                                ! Actually dimensioned (0:noSpecificChannels-1)
  END TYPE MLS_SignalsDatabase_T

CONTAINS

  ! -------------------------------------------------------------------

  ! This routine reads the signals database file and fills the data structure
  ! The routine doesn't do very much error checking as this file is assumed not
  ! to change very often.

  SUBROUTINE ReadSignalsDatabase(unit,database)

    ! Arguments and result

    INTEGER, INTENT(IN) :: unit
    TYPE (MLS_SignalsDatabase_T) :: database

    ! Local parameters ------------------------------

    INTEGER, PARAMETER :: SDB_LineLen=132
    CHARACTER (LEN=*), PARAMETER :: EOFMessage= &
         & "Unexpected EOF on signals database file"
    CHARACTER (LEN=*), PARAMETER :: AllocateMessage= &
         & "ALLOCATE failed on: "
    CHARACTER (LEN=*), PARAMETER :: DeallocateMessage= &
         & "DEALLOCATE failed on: "

    ! Some low level variables ----------------------

    CHARACTER (LEN=SDB_LineLen) :: line,first,last,rest,section
    LOGICAL :: eof
    INTEGER :: no               ! Temporary array index
    INTEGER :: signal,testSignal ! Loop counters
    INTEGER :: radiometer, band, switch, spectrometer ! More loop counters
    INTEGER :: wordLen          ! Lenght of word
    INTEGER :: hasModifier      ! Flag
    INTEGER :: status           ! From allocate

    ! These variables are intermediate arrays to allow our database to `grow'

    TYPE (SDB_SpectrometerFamilyInfo_T), DIMENSION(:), POINTER :: &
         & tempSpectrometerFamilies
    CHARACTER (LEN=SDB_LineLen), DIMENSION(:), POINTER :: &
         validSignalNames, tempValidSignalNames

    ! These strings are a breakdown of the valid signals strings

    CHARACTER (LEN=SDB_NameLen), DIMENSION(:), ALLOCATABLE :: validRadiometer, &
         & validBand, validSwitch, validSpectrometer
    CHARACTER (LEN=SDB_NameLen), DIMENSION(:), POINTER :: radiometerNames, &
         & bandNames, switchNames, spectrometerNames
    
    ! Executable code ----------------------------------------------------

    ! The first section in the signals file describes the various types
    ! of spectrometer there can be.  First, we read a line and expect it to be
    ! `spectrometers'

    CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
    IF (eof.OR.(Capitalize(TRIM(line))/="SPECTROMETERS")) &
         & CALL MLSMessage('SPECTROMETERS expected in signals database file', &
         & error=.TRUE.)

    ! Now we go through and read the family descriptions

    database%noSpectrometerFamilies=0
    SpectrometerFamilyLoop: DO
       ! Read the line
       CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
       IF (eof) CALL MLSMessage(EOFMessage,error=.TRUE.)
       IF (Capitalize(TRIM(line))=="END") EXIT SpectrometerFamilyLoop
       
       database%noSpectrometerFamilies=database%noSpectrometerFamilies+1
       no=database%noSpectrometerFamilies
       IF (no==1) THEN
          ALLOCATE(database%spectrometerFamilies(0:no-1),STAT=status)
          IF (status/=0) CALL MLSMessage(AllocateMessage // &
               & "database%spectrometerFamilies",error=.TRUE.)
       ELSE
          ALLOCATE(tempSpectrometerFamilies(0:no-1),STAT=status)
          IF (status/=0) CALL MLSMessage(AllocateMessage // &
               & "tempSpectrometerFamilies",error=.TRUE.)
          tempSpectrometerFamilies(0:no-2)=database%spectrometerFamilies
          DEALLOCATE(database%spectrometerFamilies)
          database%spectrometerFamilies=>tempSpectrometerFamilies
       ENDIF
       
       ! Now parse the line

       CALL SplitWords(line,first,rest,last,threeWay=.TRUE.,delimiter=" ")
       database%spectrometerFamilies(no-1)%name=TRIM(first)
       READ (UNIT=rest,FMT=*) database%spectrometerFamilies(no-1)%firstChannel
       READ (UNIT=last,FMT=*) database%spectrometerFamilies(no-1)%lastChannel
       database%spectrometerFamilies(no-1)%noChannels= &
            & database%spectrometerFamilies(no-1)%lastChannel- &
            & database%spectrometerFamilies(no-1)%firstChannel+1
    END DO SpectrometerFamilyLoop

    ! The next section in the file is the list of all the valid signals
    ! We'll just read these in for the moment, and parse them in the next
    ! section

    CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
    IF (eof) CALL MLSMessage(EOFMessage,error=.TRUE.)
    IF (Capitalize(TRIM(line))/="VALID SIGNALS") CALL MLSMessage( &
         & "Valid signals expected",error=.TRUE.)
    
    database%noValidSignals=0
    ValidSignalsReadingLoop: DO
       ! Read the line
       CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
       IF (eof) CALL MLSMessage(EOFMessage,error=.TRUE.)
       IF (Capitalize(TRIM(line))=="END") EXIT ValidSignalsReadingLoop
       
       database%noValidSignals=database%noValidSignals+1
       no=database%noValidSignals
       IF (no==1) THEN
          ALLOCATE(validSignalNames(0:no-1),STAT=status)
          IF (status /= 0) CALL MLSMessage(AllocateMessage // &
               & "validSignalNames",error=.TRUE.)
       ELSE
          ALLOCATE(tempValidSignalNames(0:no-1),STAT=status)
          IF (status /= 0) CALL MLSMessage(AllocateMessage // &
               & "tempValidSignalNames",error=.TRUE.)
          tempValidSignalNames(0:no-2)=validSignalNames
          DEALLOCATE(validSignalNames)
          validSignalNames=>tempValidSignalNames
       ENDIF
       validSignalNames(no-1)=line
    END DO ValidSignalsReadingLoop

    ! The remainder of the file talks about lo frequencies etc.  We'll handle
    ! all that stuff in a minute.  First we'll take apart that information we
    ! got

    ALLOCATE(validRadiometer(0:database%noValidSignals-1),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage // &
         & "validRadiometer",error=.TRUE.)
    ALLOCATE(validBand(0:database%noValidSignals-1),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage // &
         & "validBand",error=.TRUE.)
    ALLOCATE(validSwitch(0:database%noValidSignals-1),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage // &
         & "validSwitch",error=.TRUE.)
    ALLOCATE(validSpectrometer(0:database%noValidSignals-1),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage // &
         & "validSpectrometer",error=.TRUE.)
    
    DO signal=0,database%noValidSignals-1
       ! First split into radiometer,rest,spectrometer
       CALL SplitWords(validSignalNames(signal), &
            & validRadiometer(signal), rest, validSpectrometer(signal), &
            & delimiter='.',threeWay=.TRUE.)
       ! Now split rest into band,switch
       CALL SplitWords(rest,validBand(signal),validSwitch(signal), &
            & delimiter='.')
    END DO

    ! Now we find the unique ones of each of these and enter them in our
    ! database.

    CALL GetUniqueStrings(validRadiometer,radiometerNames)
    CALL GetUniqueStrings(validBand,bandNames)
    CALL GetUniqueStrings(validSwitch,switchNames)
    CALL GetUniqueStrings(validSpectrometer,spectrometerNames)

    ! Now fill up our database

    database%noRadiometers=SIZE(radiometerNames)
    database%noBands=SIZE(bandNames)
    database%noSwitches=SIZE(switchNames)

    PRINT*,"No radiometers:",database%noRadiometers
    DO radiometer=0,database%noRadiometers-1
       PRINT*,radiometerNames(radiometer)
    ENDDO

    ! Now we'll go through the radiometers and fill up the radiometerInfo

    ALLOCATE(database%radiometerInfo(0:database%noRadiometers-1),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage // &
         & "database%radiometerInfo",error=.TRUE.)
    DO radiometer=0,database%noRadiometers-1
       database%radiometerInfo(radiometer)%name=radiometerNames(radiometer)
       PRINT*,radiometerNames(radiometer)
       CALL SplitWords(radiometerNames(radiometer), &
            & first, &
            & database%radiometerInfo(radiometer)%suffix, &
            & delimiter=':')
       database%radiometerInfo(radiometer)%prefix=first
       
       ! We've filled the name, prefix and suffix.  Now parse the prefix into
       ! number and an optional modifer

       wordLen=LEN_TRIM(first)
       last=first(wordLen:wordLen)
       hasModifier=0
       IF (LLE(last,"0").AND.(LLE(last,"9"))) hasModifier=1
       PRINT *,"First: ",first
       READ (UNIT=first(2:wordLen-hasModifier),FMT=*) &
            & database%radiometerInfo(radiometer)%number
       IF (hasModifier == 1) THEN
          database%radiometerInfo(radiometer)%modifier=last
       ELSE
          database%radiometerInfo(radiometer)%modifier=""
       END IF
    END DO

    DEALLOCATE(validSignalNames)
    DEALLOCATE(radiometerNames)
    DEALLOCATE(bandNames)
    DEALLOCATE(switchNames)
    DEALLOCATE(spectrometerNames)
  END SUBROUTINE ReadSignalsDatabase


!=============================================================================
END MODULE SignalsFile
!=============================================================================

! $Log$
! Revision 1.2  1999/11/03 03:59:26  livesey
! Transfer home
!
! Revision 1.1  1999/11/03 02:53:57  livesey
! Added f90 stuff
!
 
