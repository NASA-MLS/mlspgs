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
     INTEGER :: individualChannelIndex       ! Index into array for WF4 chans
     DOUBLE PRECISION :: centerFreqIF        ! Center i.f. frequency (MHz)
  END TYPE SDB_BandInfo_T

  ! This datatype describes a spectrometer family
  TYPE SDB_SpectrometerFamilyInfo_T
     CHARACTER (LEN=SDB_NameLen) :: name ! Name of family e.g. FB25
     INTEGER :: noSpectrometers ! Number of spectrometers in this family
     INTEGER :: noChannels      ! Number of channels in family
     INTEGER :: firstChannel    ! First channel number (e.g. 1)
     INTEGER :: lastChannel     ! Last channel number (e.g. 25)
     LOGICAL :: individual      ! If set have discrete freqs. e.g. wf4 series

     ! These two arrays are the position and width of the channels wrt. the if
     ! The arrays are actually dimensioned firstChannel:lastChannel.
     ! Units are MHz

     DOUBLE PRECISION, DIMENSION(:), POINTER :: position
     DOUBLE PRECISION, DIMENSION(:), POINTER :: width
  END TYPE SDB_SpectrometerFamilyInfo_T

  ! This small datatype describes a spectrometer
  TYPE SDB_SpectrometerInfo_T
     CHARACTER (LEN=SDB_NameLen) :: name ! Name of spectrometer
     CHARACTER (LEN=SDB_NameLen) :: fullFamily ! Full name of family eg FB25
     CHARACTER (LEN=1) :: family ! Single character family id e.g. F
     INTEGER :: familyIndex     ! Index into familyInfo database
     INTEGER :: number          ! Number within family
  END TYPE SDB_SpectrometerInfo_T

  ! This datatype describes a valid radiometer/switch/band/spectrometer
  ! combination
  TYPE SDB_ValidSignal_T
     INTEGER :: radiometerIndex ! Index into array of radiometerInfos
     INTEGER :: bandIndex       ! Index into array of bandInfos
     INTEGER :: switchIndex     ! Index into array of switch names
     INTEGER :: spectrometerIndex ! Index into spectrometer array
     INTEGER :: spectrometerFamilyIndex ! Index into array of spec. fams.
     INTEGER :: spectrometerNumber ! Note count from *one* for this!
  END TYPE SDB_ValidSignal_T

  ! This datatype is an amalgam of the above and is the database that is filled
  TYPE MLS_SignalsDatabase_T
     INTEGER :: noRadiometers   ! Including redundant etc.
     INTEGER :: noBands         ! Accross the whole instrument
     INTEGER :: noSwitches      ! No. vald S0, S1 etc. fields
     INTEGER :: noSpectrometers ! No. spectrometers in whole instrument
     INTEGER :: noSpectrometerFamilies
     INTEGER :: noValidSignals    ! Number of valid Signal combinations
     INTEGER :: noIndividualChannels ! Total no. channels in WF4 type specs.
     
     TYPE (SDB_RadiometerInfo_T), DIMENSION(:), POINTER :: radiometerInfo
                                ! Actually dimensioned (0:noRadiometers-1)
     TYPE (SDB_BandInfo_T), DIMENSION(:), POINTER :: bandInfo
                                ! Actually dimensioned (0:noBands-1)
     CHARACTER (LEN=SDB_NameLen), DIMENSION(:), POINTER :: switches
                                ! Actually dimensioned (0:noSwitches-1)
     TYPE (SDB_SpectrometerInfo_T), DIMENSION(:), POINTER :: spectrometerInfo
     TYPE (SDB_SpectrometerFamilyInfo_T), DIMENSION(:), POINTER :: &
          & spectrometerFamilies
             ! Actually dimensioned (0:noSpectrometerFamilies-1)
     TYPE (SDB_ValidSignal_T), DIMENSION(:), POINTER :: validSignals
                                ! Valid Signal combinations
                                ! Actually dimensioned (0:noValidSignals-1)

     ! These two are for `individual' channels, they describe their position and
     ! width in i.f. space in MHz.
     DOUBLE PRECISION, DIMENSION(:), POINTER :: individualPosition, &
          & individualWidth
                            ! Actually dimensioned (0:noIndividualChannels-1)
  END TYPE MLS_SignalsDatabase_T

CONTAINS

  ! -------------------------------------------------------------------

  ! This routine reads the signals database file and fills the data structure
  ! The routine doesn't do very much error checking as this file is assumed not
  ! to change very often.

  SUBROUTINE ReadSignalsDatabase(unit,database)

    ! Arguments and result

    INTEGER, INTENT(IN) :: unit
    TYPE (MLS_SignalsDatabase_T), INTENT(OUT) :: database

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
    INTEGER :: spectrometerFamily ! Another loop counter
    INTEGER :: wordLen          ! Length of word
    INTEGER :: hasModifier      ! Flag
    INTEGER :: status           ! From allocate
    INTEGER :: index            ! General array index

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

    CHARACTER (LEN=1), DIMENSION(:), ALLOCATABLE :: &
         & spectrometerFamilyChars

    INTEGER :: evenNo           ! For `even' channels
    DOUBLE PRECISION :: evenStart,evenSpacing,evenWidth ! For `even' channels
    
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

    ALLOCATE(radiometerNames(0:database%noValidSignals-1),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage // &
         & "radiometerNames",error=.TRUE.)
    ALLOCATE(bandNames(0:database%noValidSignals-1),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage // &
         & "bandNames",error=.TRUE.)
    ALLOCATE(switchNames(0:database%noValidSignals-1),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage // &
         & "switchNames",error=.TRUE.)
    ALLOCATE(spectrometerNames(0:database%noValidSignals-1),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage // &
         & "spectrometerNames",error=.TRUE.)
    
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

    CALL GetUniqueStrings(validRadiometer,radiometerNames, &
         & database%noRadiometers)
    CALL GetUniqueStrings(validBand,bandNames,database%noBands)
    CALL GetUniqueStrings(validSwitch,switchNames,database%noSwitches)
    CALL GetUniqueStrings(validSpectrometer,spectrometerNames, &
         & database%noSpectrometers)

    ! Now fill up our database element by element

    ! We'll go through the radiometers and fill up the radiometerInfo

    ALLOCATE(database%radiometerInfo(0:database%noRadiometers-1),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage // &
         & "database%radiometerInfo",error=.TRUE.)

    DO radiometer=0,database%noRadiometers-1
       database%radiometerInfo(radiometer)%name=radiometerNames(radiometer)
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
       IF (LLT(last,"0").OR.(LGT(last,"9"))) hasModifier=1
       READ (UNIT=first(2:wordLen-hasModifier),FMT=*) &
            & database%radiometerInfo(radiometer)%number
       IF (hasModifier == 1) THEN
          database%radiometerInfo(radiometer)%modifier=last
       ELSE
          database%radiometerInfo(radiometer)%modifier=""
       END IF
    END DO

    ! Now, we similarly go through the bands and fill up that info

    ALLOCATE(database%bandInfo(0:database%noBands-1),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage // &
         & "database%bandInfo",error=.TRUE.)

    DO band=0,database%noBands-1
       database%bandInfo(band)%name=bandNames(band)
       CALL SplitWords(bandNames(band), &
            & first, &
            & database%bandInfo(band)%suffix, &
            & delimiter=':')

       ! Now we parse the prefix section.  This is a B followed by a number
       ! then a spectrometer type
       
       wordLen=LEN_TRIM(first)
       database%bandInfo(band)%spectrometerFamily=first(wordLen:wordLen)
       READ (UNIT=first(2:wordLen-1),FMT=*) &
            database%bandInfo(band)%number
    END DO ! We'll sort out the spectrometer family index later

    ! Sorting out the switch settings is easy.

    ALLOCATE (database%switches(0:database%noSwitches-1))
    database%switches=switchNames

    ! Finally we do the spectrometers

    ALLOCATE (database%spectrometerInfo(0:database%noSpectrometers-1), &
         &  STAT=status)
    IF (status /= 0) CALL MLSMessage( &
         & AllocateMessage//"database%spectrometerInfo",error=.TRUE.)
    DO spectrometer=0,database%noSpectrometers-1
       database%spectrometerInfo(spectrometer)%name= &
            & spectrometerNames(spectrometer)
       CALL SplitWords(spectrometerNames(spectrometer), &
            & database%spectrometerInfo(spectrometer)%fullFamily, &
            rest,delimiter="-")
       READ (UNIT=rest,FMT=*) &
            database%spectrometerInfo(spectrometer)%number
       database%spectrometerInfo(spectrometer)%family= &
            spectrometerNames(spectrometer)(1:1)       
    END DO

    ! Now as a final tidy up, we go through and fill in the family index
    ! variables in both the spectrometerInfo and bandInfo database entries.

    ALLOCATE (spectrometerFamilyChars(0:database%noSpectrometerFamilies-1))

    DO spectrometerFamily=0,database%noSpectrometerFamilies-1
       spectrometerFamilyChars(spectrometerFamily)= &
            & database%spectrometerFamilies(spectrometerFamily)%name(1:1)
    END DO

    DO spectrometer=0,database%noSpectrometers-1
       database%spectrometerInfo(spectrometer)%familyIndex= &
            & LinearSearchStringArray(database%spectrometerFamilies%name, &
            & database%spectrometerInfo(spectrometer)%fullFamily)-1          
    END DO

    DO band=0,database%noBands-1
       database%bandInfo(band)%spectrometerFamilyIndex= &
            & LinearSearchStringArray(spectrometerFamilyChars, &
            & database%bandInfo(band)%spectrometerFamily)
    END DO

    ! That's basically the valid signals database filled up, no we go on and
    ! read the rest of the file which contains frequency information etc.

    ! The next section descusses the radiometer frequencies
    
    CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
    IF (eof.OR.(Capitalize(TRIM(line))/="RADIOMETER FREQUENCIES")) &
         & CALL MLSMessage("RADIOMETER FREQUENCIES expected in database", &
         & error=.TRUE.)

    ! Now this section is just a list of names followed by lo's

    database%radiometerInfo%lo=0.0D0

    RadiometerFrequencyLoop: DO
       CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
       IF (eof) CALL MLSMessage(EOFMessage,error=.TRUE.)
       IF (Capitalize(TRIM(line))=="END") EXIT RadiometerFrequencyLoop
       CALL SplitWords(line,first,rest,delimiter=' ')
       radiometer=LinearSearchStringArray(database%radiometerInfo%name, &
            & first, caseInsensitive=.TRUE.)-1
       ! Note the -1 due to parameter paassing.
       IF (radiometer == -1) CALL MLSMessage("No such radiometer: "//first, &
            & error=.TRUE.)
       READ (unit=rest,FMT=*) database%radiometerInfo(radiometer)%lo
    END DO RadiometerFrequencyLoop

    IF (MINVAL(database%radiometerInfo%lo) <= 0.0D0) CALL MLSMessage( &
         & "Not all radiometer LOs assigned",error=.TRUE.)

    ! Now the next section is a set of similar information for the bands

    CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
    IF (eof.OR.(Capitalize(TRIM(line))/="BAND FREQUENCIES")) &
         & CALL MLSMessage("BAND FREQUENCIES expected in database", &
         & error=.TRUE.)

    ! Now this section is just a list of names followed by center frequencies

    database%bandInfo%centerFreqIF=0.0D0

    BandFrequencyLoop: DO
       CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
       IF (eof) CALL MLSMessage(EOFMessage,error=.TRUE.)
       IF (Capitalize(TRIM(line))=="END") EXIT BandFrequencyLoop
       CALL SplitWords(line,first,rest,delimiter=" ")
       band=LinearSearchStringArray(database%bandInfo%name, &
            & first, caseInsensitive=.TRUE.)-1
       ! Note -1 due to parameter passing
       IF (band == -1) CALL MLSMessage("No such band: "//first, &
            & error=.TRUE.)
       READ (unit=rest,FMT=*) database%bandInfo(band)%centerFreqIF
    END DO BandFrequencyLoop

    ! Don't check here for all allocated as the WF type filters are listed
    ! seperately.

    ! Now we have the section giving the frequencies and widths of all the
    ! spectrometer channels.

    CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
    IF (eof.OR.(Capitalize(TRIM(line))/="SPECTROMETER FREQUENCIES")) &
         & CALL MLSMessage("SPECTROMETER FREQUENCIES expected in database", &
         & error=.TRUE.)

    SpectrometerFrequencyLoop: DO
       CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
       IF (eof) CALL MLSMessage(EOFMessage,error=.TRUE.)
       IF (Capitalize(TRIM(line))=="END") EXIT SpectrometerFrequencyLoop

       ! Each spectrometer family has a line describing it with a set of
       ! instructions as to how to assign the channels

       CALL SplitWords(line,first,rest,delimiter=" ")
       index=LinearSearchStringArray( &
            & database%spectrometerFamilies%name,first, &
            & caseInsensitive=.TRUE.)-1
       ! Note the -1 due to parameter passing
       IF (index == -1) CALL MLSMessage( &
            & "No such spectrometer family: "//first,error=.TRUE.)

       ! Here we'll allocate the arrays for size and posiiton.  Note for the WF
       ! filters this is unneccessary, and we'll deallocate them later.
       ! But to make the code more readable I'm allocating them first.

       ALLOCATE (database%spectrometerFamilies(index)%position(&
        & database%spectrometerFamilies(index)%firstChannel:&
        & database%spectrometerFamilies(index)%lastChannel),&
        & STAT=status)
       IF (status /= 0) CALL MLSMessage(AllocateMessage//first//" position")

       ALLOCATE (database%spectrometerFamilies(index)%width(&
        & database%spectrometerFamilies(index)%firstChannel:&
        & database%spectrometerFamilies(index)%lastChannel),&
        & STAT=status)
       IF (status /= 0) CALL MLSMessage(AllocateMessage//first//" width")

       database%spectrometerFamilies(index)%individual=.FALSE.

       SELECT CASE (Capitalize(TRIM(rest)))

          CASE ('LIST') ! Positions and widths given on next lines
             READ (UNIT=unit, FMT=*) &
                  & database%spectrometerFamilies(index)%position
             READ (UNIT=unit, FMT=*) &
                  & database%spectrometerFamilies(index)%width

          CASE ('EVEN') ! No chans, start Freq, spacing and widths
             READ (UNIT=unit, FMT=*) evenNo,evenStart,evenSpacing,evenWidth
             IF (evenNo /= database%spectrometerFamilies(index)%noChannels) &
                  & CALL MLSMessage("Wrong number of channels for "//first, &
                  & error=.TRUE.)
             ! Loop and fill information
             database%spectrometerFamilies(index)%width=evenWidth
             DO evenNo=database%spectrometerFamilies(index)%firstChannel, &
                   & database%spectrometerFamilies(index)%lastChannel
                database%spectrometerFamilies(index)%position= &
                     & evenStart+evenNo*evenSpacing
             END DO

             CASE ('INDIVIDUAL')
                database%spectrometerFamilies(index)%individual=.FALSE.
                DEALLOCATE(database%spectrometerFamilies(index)%position)
                DEALLOCATE(database%spectrometerFamilies(index)%width)

             CASE DEFAULT 
                CALL MLSMessage("Unrecognized spectrometer family: "//first, &
                     & error=.TRUE.)
       END SELECT
    END DO SpectrometerFrequencyLoop

    ! The final section of the file is the section dealing with specific
    ! channels. Here, having built up most of the database, we can find the
    ! relevant data using our standard parse routines.

    ! Now we have read all the gory information from our database
    ! we need to fill up our valid signals database

    ALLOCATE (database%validSignals(0:database%noValidSignals-1),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage//"database%validSignals", &
         & error=.TRUE.)

    DO signal=0,database%noValidSignals-1
       ! Note here we have to -1 in each case because of the parameter passing
       ! issue.

       database%validSignals(signal)%radiometerIndex=LinearSearchStringArray( &
            & database%radiometerInfo%name,validRadiometer(signal))-1

       database%validSignals(signal)%bandIndex=LinearSearchStringArray( &
            & database%bandInfo%name,validBand(signal))-1

       database%validSignals(signal)%switchIndex=LinearSearchStringArray( &
            & database%switches,validSwitch(signal))

       ! For the spectrometers we need to do more
       index=LinearSearchStringArray(database%spectrometerInfo%name, &
            & validSpectrometer(signal))-1

       database%validSignals(signal)%spectrometerIndex=index
       database%validSignals(signal)%spectrometerFamilyIndex= &
            & database%spectrometerInfo(index)%familyIndex
       database%validSignals(signal)%spectrometerNumber= &
            & database%spectrometerInfo(index)%number
    END DO

    ! Now we tidy up our arrays and exit

    DEALLOCATE(validSignalNames)
    DEALLOCATE(radiometerNames)
    DEALLOCATE(bandNames)
    DEALLOCATE(switchNames)
    DEALLOCATE(spectrometerNames)
    DEALLOCATE(spectrometerFamilyChars)
  END SUBROUTINE ReadSignalsDatabase

  ! -------------------------------------------------------------------

  ! This routine deallocates all the information dealing with the signals
  ! database

  SUBROUTINE DestroySignalsDatabase(database)

    ! Argument

    TYPE (MLS_SignalsDatabase_T), INTENT(INOUT) :: database
    
    ! Local variable

    INTEGER :: i

    ! Executable code

    database%noRadiometers=0
    database%noBands=0
    database%noSwitches=0
    database%noSpectrometers=0
    database%noSpectrometerFamilies=0
    database%noValidSignals=0
    database%noIndividualChannels=0

    DEALLOCATE (database%radiometerInfo)
    DEALLOCATE (database%bandInfo)
    DEALLOCATE (database%switches)
    DEALLOCATE (database%spectrometerInfo)

    DO i=0,database%noSpectrometerFamilies-1
       DEALLOCATE (database%spectrometerFamilies(i)%position)
       DEALLOCATE (database%spectrometerFamilies(i)%width)
    END DO
    DEALLOCATE (database%spectrometerFamilies)

    DEALLOCATE (database%validSignals)
    DEALLOCATE (database%individualPosition)
    DEALLOCATE (database%individualWidth)

  END SUBROUTINE DestroySignalsDatabase

!=============================================================================
END MODULE SignalsFile
!=============================================================================

! $Log$
! Revision 1.4  1999/11/04 00:56:34  livesey
! Transfer home
!
! Revision 1.3  1999/11/04 00:06:43  livesey
! Removed the irrelevant old history from the other day.
!
! Revision 1.2  1999/11/03 23:56:26  livesey
! A trivial change to test out the cvs watch stuff.
!
! Revision 1.1  1999/11/03 23:53:58  livesey
! Added mlsstrings and signalsfile
 
