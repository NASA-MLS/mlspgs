! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7407 is acknowledged.

!=============================================================================
MODULE MLSSignalNomenclature    ! Dealing with MLS rad.band etc. specifiers
!=============================================================================

  USE MLSCommon
  USE MLSStrings
  USE MLSMessageModule
  USE Intrinsic, ONLY: L_None

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: MLSSignal_T, ParseMLSSignalRequest, DestroyMLSSignalsInfo,&
       & ConcatenateMLSSignalsInfo, UnionMLSSignalsInfo, &
       & IntersectionMLSSignalsInfo, GetFullMLSSignalName, &
       & ReadSignalsDatabase, DestroySignalsDatabase, &
       & GetMLSRadiometerNames, GetMLSBandNames

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = & 
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  ! The same as 2.1 before being removed (see log below)

  ! This module contains all the functionality required for dealing with the
  ! standard MLS signal nomenclature.  The instrument configuration is read
  ! from a file (e.g. emls-signals.dat). The user can request information about
  ! a signal or set of signals by name.

  ! This datatype describes a valid radiometer/switch/band/spectrometer
  ! combination.  It is the one most used by calling code.  For this reason I
  ! have put it first.  Some of the information given is duplicated in later
  ! structures but they are typically only internally used by this module.

  TYPE MLSSignal_T

     ! First we have a set of indices into arrays.  These are used a lot in the
     ! module.  But the calling code may also find them useful for issues such
     ! as: Does signal A come from the same radiometer as signal B?

     INTEGER :: signalDatabaseIndex ! Index into rad.band.s.spec combination
     INTEGER :: instrumentModule ! Module in the instrument (EMLS:1=GHz, 2=THz)
     INTEGER :: radiometerIndex ! Index into array of radiometerInfos
     INTEGER :: bandIndex       ! Index into array of bandInfos
     INTEGER :: switchIndex     ! Index into array of switch names
     INTEGER :: spectrometerIndex ! Index into spectrometer array
     INTEGER :: upperLower      ! -1=Lower, 1=Upper, 0=Folded
     LOGICAL :: notCopy         ! If set POINTER arrays below are not copies

     ! For an explanation of the notCopy parameter see the
     ! ParseMLSSignalRequest routine.

     ! Now we have detailed information on the radiometer

     CHARACTER (LEN=NameLen) :: radiometerName   ! e.g. R1A:118
     CHARACTER (LEN=NameLen) :: radiometerPrefix ! e.g. R1A
     CHARACTER (LEN=NameLen) :: radiometerSuffix ! e.g. 118
     INTEGER :: radiometerNumber                    ! e.g. 1
     CHARACTER (LEN=1) :: radiometerModifier        ! e.g. A
     REAL(r8) :: lo                         ! MHz

     ! Now we have detailed information on the band

     CHARACTER (LEN=NameLen) :: bandName    ! e.g. B1F:PT
     CHARACTER (LEN=NameLen) :: bandSuffix  ! e.g. PT
     CHARACTER (LEN=1) :: spectrometerFamily   ! e.g. F
     REAL(r8) :: bandCenterFreqIF      ! MHz

     ! Now we have information on the switch

     CHARACTER (LEN=NameLen) :: switch

     ! Now information on the spectrometers

     CHARACTER (LEN=NameLen) :: spectrometerName
     CHARACTER (LEN=NameLen) :: fullSpectrometerFamily ! e.g. FB25
     INTEGER :: spectrometerFamilyIndex ! Index into array of spec. fams.
     INTEGER :: spectrometerNumber ! Note count from one for this.

     ! Now information the channels in the whole band

     INTEGER :: firstChannelInBand
     INTEGER :: lastChannelInBand
     INTEGER :: noChannelsInBand

     ! Now this particular signal may be a subset of all the channels so
     ! detail what we have.

     INTEGER :: noChannelsIncluded
     LOGICAL, DIMENSION(:), POINTER :: channelIncluded
     REAL(r8), DIMENSION(:), POINTER :: channelPosition ! i.f. space
     REAL(r8), DIMENSION(:), POINTER :: channelWidth ! i.f. space

  END TYPE MLSSignal_T

  ! ----------------------------------------------------------------------

  ! The remaining datatypes are somewhat more private.

  ! This datatype describes a radiometer
  TYPE SDBRadiometerInfo_T
     CHARACTER (LEN=NameLen) :: name   ! e.g. R1A:118
     CHARACTER (LEN=NameLen) :: prefix ! e.g. R1A
     CHARACTER (LEN=NameLen) :: suffix ! e.g. 118
     INTEGER :: number                     ! e.g. 1
     CHARACTER(LEN=1) :: modifier          ! e.g. A/B or H/V for R5 (emls)
     REAL(r8) :: lo                ! Local oscillator /MHz
     INTEGER :: instrumentModule ! Module in instrument GHz/THz 
  END TYPE SDBRadiometerInfo_T

  ! This datatype describes a band
  TYPE SDBBandInfo_T
     CHARACTER (LEN=NameLen) :: name      ! e.g. B1F:PT
     CHARACTER (LEN=NameLen) :: suffix    ! e.g. PT
     INTEGER :: number                       ! e.g. 1
     CHARACTER (LEN=1) :: spectrometerFamily ! e.g. F
     INTEGER :: spectrometerFamilyIndex      ! Index into array of next type
     REAL(r8) :: centerFreqIF        ! Center i.f. frequency (MHz)
  END TYPE SDBBandInfo_T

  ! This datatype describes a spectrometer family
  TYPE SDBSpectrometerFamilyInfo_T
     CHARACTER (LEN=NameLen) :: name ! Name of family e.g. FB25
     INTEGER :: noSpectrometers ! Number of spectrometers in this family
     INTEGER :: noChannels      ! Number of channels in family
     INTEGER :: firstChannel    ! First channel number (e.g. 1)
     INTEGER :: lastChannel     ! Last channel number (e.g. 25)
     LOGICAL :: individual      ! If set have discrete freqs. e.g. wf4 series

     ! These two arrays are the position and width of the channels wrt. the if
     ! The arrays are actually dimensioned firstChannel:lastChannel.
     ! Units are MHz

     REAL(r8), DIMENSION(:), POINTER :: position
     REAL(r8), DIMENSION(:), POINTER :: width
  END TYPE SDBSpectrometerFamilyInfo_T

  ! This small datatype describes a spectrometer
  TYPE SDBSpectrometerInfo_T
     CHARACTER (LEN=NameLen) :: name ! Name of spectrometer
     CHARACTER (LEN=NameLen) :: fullFamily ! Full name of family eg FB25
     CHARACTER (LEN=1) :: family ! Single character family id e.g. F
     INTEGER :: familyIndex     ! Index into familyInfo database
     INTEGER :: number          ! Number within family
  END TYPE SDBSpectrometerInfo_T

  ! This datatype is an amalgam of the above and is the database that is filled
  TYPE MLSSignalsDatabase_T
     INTEGER :: noRadiometers   ! Including redundant etc.
     INTEGER :: noBands         ! Accross the whole instrument
     INTEGER :: noSwitches      ! No. vald S0, S1 etc. fields
     INTEGER :: noSpectrometers ! No. spectrometers in whole instrument
     INTEGER :: noSpectrometerFamilies
     INTEGER :: noValidSignals    ! Number of valid Signal combinations

     TYPE (SDBRadiometerInfo_T), DIMENSION(:), POINTER :: radiometerInfo
          ! Actually dimensioned (noRadiometers)
     TYPE (SDBBandInfo_T), DIMENSION(:), POINTER :: bandInfo
          ! Actually dimensioned (noBands)
     CHARACTER (LEN=NameLen), DIMENSION(:), POINTER :: switches
          ! Actually dimensioned (noSwitches)
     TYPE (SDBSpectrometerInfo_T), DIMENSION(:), POINTER :: spectrometerInfo
     TYPE (SDBSpectrometerFamilyInfo_T), DIMENSION(:), POINTER :: &
          & spectrometerFamilyInfo
          ! Actually dimensioned (noSpectrometerFamilies)
     TYPE (MLSSignal_T), DIMENSION(:), POINTER :: validSignals
          ! Valid Signal combinations
          ! Actually dimensioned (noValidSignals)
  END TYPE MLSSignalsDatabase_T

  ! Local, private variables

  TYPE (MLSSignalsDatabase_T) :: database

CONTAINS

  ! ===================================================================

  ! The first routine parses a request for a particular radiometer
  ! and returns an array of flags indicating which radiometers match

  SUBROUTINE ParseRadiometerRequest (request,matches)

    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN) :: request

    LOGICAL, DIMENSION(:) :: matches ! Result (database%noRadiometers)

    ! Local variables
    CHARACTER (LEN=LEN(request)) :: prefix,suffix
    CHARACTER (LEN=1) :: modifierRequested
    INTEGER :: numberRequested,hasModifier,prefixLen,radiometer

    ! Executable code

    CALL SplitWords(request,prefix,suffix,delimiter=':')
    prefix=Capitalize(prefix)
    suffix=Capitalize(suffix)

    IF (prefix(1:1)/="R") CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "R expected in radiometer specifier")

    ! Now parse the rest of the prefix.
    ! The format of the field is R<number><modifier>, where modifier may be
    ! omitted or set to *. Also R* is valid, or just R!

    prefixLen=LEN_TRIM(prefix)
    numberRequested=0
    modifierRequested="*"
    hasModifier=1
    
    ! First we'll look for the modifier, look at the last characater.  If it's
    ! numerical set the modifier to "*"

    IF (prefixLen>1) THEN
       modifierRequested=prefix(prefixLen:prefixLen)
       IF ((LGE(modifierRequested,"0")).AND.(LLE(modifierRequested,"9"))) THEN
          hasModifier=0
          modifierRequested="*"
       ENDIF

       ! Now we look at the numeric field if there is one.

       IF (prefixLen-hasModifier > 1) &
            & READ (UNIT=prefix(2:prefixLen-hasModifier),FMT=*) numberRequested
    ENDIF

    ! Now we have the requsted number (or 0 if dont care) and requested
    ! modifier (or * if dont care)

    IF (SIZE(matches) /= database%noRadiometers) CALL MLSMessage( &
         MLSMSG_Error,ModuleName,"Result is wrong size")

    matches= &
         & ( (numberRequested == database%radiometerInfo%number) .OR. &
         &   (numberRequested == 0) ) .AND. &
         & ( (modifierRequested == database%radiometerInfo%modifier) .OR. &
         &   (modifierRequested == "*") )

    ! To match the suffix we have to be a little more old fasioned because we
    ! want to use trim.

    IF (TRIM(suffix) /= "") THEN
       DO radiometer=1,database%noRadiometers
          matches(radiometer) = matches(radiometer).AND.&
            & (TRIM(suffix) ==TRIM(database%radiometerInfo(radiometer)%suffix))
       ENDDO
    ENDIF

  END SUBROUTINE ParseRadiometerRequest

  ! -------------------------------------------------------------------

  ! This function is similar to the above one, except that more can be omitted
  ! (e.g. the spectrometer type etc.) Also the user can specify a request for
  ! upper/lower sideband.  This is reflected in the upper/lower value
  ! -1=lower, 1=upper, 0=folded.

  SUBROUTINE ParseBandRequest (request,matches,upperLower)

    ! Dummy arguements
    CHARACTER (LEN=*), INTENT(IN) :: request
    LOGICAL, DIMENSION(:) :: matches ! Result (database%noBands)
    INTEGER, INTENT(OUT), OPTIONAL :: upperLower

    ! Local variables
    CHARACTER (LEN=LEN(request)) :: prefix,suffix
    CHARACTER (LEN=1) :: thisChar, spectrometerFamilyRequested
    INTEGER :: upperLowerRequested
    INTEGER :: pos,bandNumberRequested,prefixLen,band
    LOGICAL :: bandStarRequested

    ! Executable code

    CALL SplitWords(request,prefix,suffix,delimiter=":")
    prefix=Capitalize(prefix)
    suffix=Capitalize(suffix)
    thisChar=prefix(1:1)
    IF (thisChar /= "B") CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "B expected in band specifier")

    ! Start with a clean slate, the user is not fussy

    upperLowerRequested=0
    spectrometerFamilyRequested="*"
    bandNumberRequested=0
    bandStarRequested=.FALSE.

    prefixLen=LEN_TRIM(prefix)
    IF (prefixLen>1) THEN

       ! OK so it's not just B: something, it's more complicated

       ! Look at the second character, it's either a number or a *
       thisChar=prefix(2:2)
       IF (thisChar == "*") bandStarRequested=.TRUE.

       ! Now loop in from the back end to take the prefix apart
       pos=prefixLen
       ParseBandRequestParse: DO
          thisChar=prefix(pos:pos)
          IF (LGE(thisChar,"0").AND.(LLE(thisChar,"9"))) &
               & EXIT ParseBandRequestParse
          SELECT CASE (thisChar)
          CASE ("*")
             IF (bandStarRequested.AND.(pos==2)) EXIT ParseBandRequestParse
             ! Otherwise it's a spectrometer family request, which is already
             ! set to * so there's nothing to do.
          CASE ("U")
             upperLowerRequested=1
          CASE ("L")
             upperLowerRequested=-1
          CASE DEFAULT
             spectrometerFamilyRequested=thisChar
          END SELECT
          pos=pos-1
          IF (pos == 1) EXIT ParseBandRequestParse
       END DO ParseBandRequestParse
       IF (pos==1) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & "Bad band request: "//TRIM(request))
       
       ! Now read the band number from prefix
       
       IF (.NOT. bandStarRequested) READ (UNIT=prefix(2:pos),FMT=*) &
            & bandNumberRequested

    ENDIF                       ! Not just B:...

    ! Now we have bandNumberRequested or 0 for don't care,
    ! UpperLower selected and spectrometerFamilyRequest
    ! Find all the bands that match our criteria.

    IF (SIZE(matches) /= database%noBands) CALL MLSMessage(&
         & MLSMSG_Error,ModuleName,"Result wrong size")

    matches= &
         & ( (bandNumberRequested == database%bandInfo%number) .OR. &
         &   (bandNumberRequested == 0) ) .AND. &
         & ( (spectrometerFamilyRequested == &
         &       database%bandInfo%spectrometerFamily) .OR. &
         &   (spectrometerFamilyRequested == "*") )

    ! For the suffix we have to be a little more old fasioned because we can't
    ! do a `ragged' trim

    IF (TRIM(suffix) /= "") THEN
       DO band=1,database%noBands
          matches(band)=matches(band) .AND. &
               & (TRIM(suffix) == TRIM(database%bandInfo(band)%suffix))
       ENDDO
    ENDIF

    IF (PRESENT(upperLower)) upperLower=upperLowerRequested
  END SUBROUTINE ParseBandRequest

  ! -------------------------------------------------------------------

  ! This routine parses a switch request

  SUBROUTINE ParseSwitchRequest(request,matches)

    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN) :: request
    LOGICAL, DIMENSION(:) :: matches ! (database%noSwitches)

    ! Local variables
    INTEGER :: switch
    CHARACTER (LEN=LEN(request)) :: capRequest

    ! Executable code

    capRequest=Capitalize(request)
    IF (capRequest(1:1) /= "S") CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "S expected for switch spec")
    IF (SIZE(matches) /= database%noSwitches) CALL MLSMessage( &
         MLSMSG_Error,ModuleName,"Result wrong size")

    DO switch=1,database%noSwitches
       matches(switch)=(capRequest==database%switches(switch))
    ENDDO
  END SUBROUTINE ParseSwitchRequest

  ! -------------------------------------------------------------------

  ! This routine parses a spectrometer request. There's no need for wildcards
  ! here that make any sense.

  SUBROUTINE ParseSpectrometerRequest(request,matches)

    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN) :: request
    LOGICAL, DIMENSION(:) :: matches ! (database%noSpectrometers)

    ! Local variables
    CHARACTER (LEN=LEN(request)) :: capRequest
    INTEGER :: spectrometer

    ! Executable code

    IF (SIZE(matches) /= database%noSpectrometers) CALL MLSMessage( &
         & MLSMSG_Error,ModuleName,"Result wrong size")

    capRequest=Capitalize(request)
    DO spectrometer=1,database%noSpectrometers
       matches(spectrometer) = (TRIM(capRequest) == &
            & TRIM(database%spectrometerInfo(spectrometer)%name))
    ENDDO

  END SUBROUTINE ParseSpectrometerRequest
  
  ! -------------------------------------------------------------------

  ! This function parses a channel request.  The form is 
  ! C[<spec>+<spec>+<spec>] etc. where <spec> is n or m:n

  SUBROUTINE ParseChannelRequest ( request,firstChannel,lastChannel, &
    & channelIncluded )

    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN) :: request
    INTEGER, INTENT(IN) :: firstChannel, lastChannel
    LOGICAL, DIMENSION(:), POINTER :: channelIncluded

    ! Local variables
    CHARACTER (LEN=LEN(request)) :: field, remainder, newRemainder, &
      & firstStr, lastStr
    INTEGER :: highestChannel, firstInRange, lastInRange, channel
    INTEGER :: status, requestLen
    CHARACTER (LEN=1) :: keyChar ! Needed due to absoft bug

    ! Executable code

    highestChannel=firstChannel-1
    keyChar=request(1:1)
    IF ((keyChar/="C").AND.(keyChar/="c")) CALL MLSMessage(MLSMSG_Error, &
         & ModuleName,"C expected in channel specifier")
    keyChar=request(2:2)
    IF (keyChar/="[") CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "[ expected in channel specifier")
    requestLen=LEN_TRIM(request)
    keyChar=request(requestLen:requestLen)
    IF (keyChar/="]") CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "] expected in channel specifier")

    remainder=Capitalize(request(3:requestLen-1))
    ALLOCATE(channelIncluded(firstChannel:lastChannel),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Allocate failed for channelIncluded")
    channelIncluded=.FALSE.

    ParseChannelParseLoop: DO
       CALL SplitWords(remainder,field,newRemainder,delimiter="+")
       remainder=newRemainder
       IF (TRIM(field)=="") EXIT ParseChannelParseLoop

       ! See if this is a range or just a single channel
       IF (INDEX(field,":")/=0) THEN
          ! It's a range
          CALL SplitWords(field,firstStr,lastStr,delimiter=":")

          IF (TRIM(firstStr)=="") THEN 
             firstInRange=firstChannel
          ELSE
             READ (UNIT=firstStr,FMT=*) firstInRange
          ENDIF

          IF (TRIM(lastStr)=="") THEN
             lastInRange=lastChannel
          ELSE
             READ (UNIT=lastStr,FMT=*) lastInRange
          ENDIF
       ELSE
          ! It's just a single channel
          READ (UNIT=field,FMT=*) firstInRange
          lastInRange=firstInRange
       ENDIF

       IF (firstInRange<=highestChannel) CALL MLSMessage(MLSMSG_Error, &
            & ModuleName, "Channel specifier out of order")
       DO channel=firstInRange,lastInRange
          channelIncluded(channel)=.TRUE.
       ENDDO
       highestChannel=lastInRange
    END DO ParseChannelParseLoop

  END SUBROUTINE ParseChannelRequest

  ! -------------------------------------------------------------------

  ! This routine makes sure that the channel information in a signal is a copy
  ! of the original, and doesn't merely point to it.

  SUBROUTINE TurnMLSChannelInfoIntoCopy(signals)

    ! Dummy argument
    TYPE (MLSSignal_T), DIMENSION(:), INTENT(INOUT) :: signals

    ! Local variables
    INTEGER :: status, signal
    REAL(r8), DIMENSION(:), POINTER :: tempPosition, tempWidth
    LOGICAL, DIMENSION(:), POINTER :: tempIncluded

    DO signal=1,SIZE(signals)
       tempPosition=>signals(signal)%channelPosition
       tempWidth=>signals(signal)%channelWidth
       tempIncluded=>signals(signal)%channelIncluded
       
       ALLOCATE(signals(signal)%channelPosition( &
            & signals(signal)%firstChannelInBand: &
            & signals(signal)%lastChannelInBand),STAT=status)
       IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"channelPosition")
       ALLOCATE(signals(signal)%channelWidth( &
            & signals(signal)%firstChannelInBand: &
            & signals(signal)%lastChannelInBand),STAT=status)
       IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"channelWidth")
       ALLOCATE(signals(signal)%channelIncluded( &
            & signals(signal)%firstChannelInBand: &
            & signals(signal)%lastChannelInBand),STAT=status)
       IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"channelIncluded")

       signals(signal)%channelPosition=tempPosition
       signals(signal)%channelWidth=tempWidth
       signals(signal)%channelIncluded=tempIncluded
    END DO
  END SUBROUTINE TurnMLSChannelInfoIntoCopy

  ! -------------------------------------------------------------------

  ! This subroutine takes a full request and returns an index or set of
  ! indices into the valid signals data structure, along with an indication of
  ! upper/lower sideband if relevant.  Also returned is a list of channels.

  SUBROUTINE ParseMLSSignalRequest(request, &
       & signals,noCopy)

    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN) :: request
    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signals
    LOGICAL, OPTIONAL, INTENT(IN) :: noCopy ! Don't copy data, point to it

    ! Local paramters
    CHARACTER (LEN=*), PARAMETER :: OutOfOrder="Signal designation out of order"

    ! Local variables

    INTEGER :: upperLower
    INTEGER :: mostAdvancedField
    INTEGER :: radiometerIndex,bandIndex,switchIndex,spectrometerIndex
    CHARACTER (LEN=LEN(request)) :: remainder,newRemainder,field, &
         & channelString
    LOGICAL, DIMENSION(:), ALLOCATABLE :: &
         & radiometerMatches, bandMatches, switchMatches, &
         & spectrometerMatches, allMatch
    INTEGER :: noMatches,status,signal

    LOGICAL, DIMENSION(:), POINTER :: channelIncluded
    LOGICAL :: uniqueSpectrometerFamilyRequest
    CHARACTER (LEN=1) :: spectrometerFamily
    LOGICAL :: useNoCopy
    CHARACTER (LEN=1) :: keyChar ! Needed due to absoft bug

    ! Executable code

    ALLOCATE (radiometerMatches(database%noRadiometers),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"radiometerMatches")
    ALLOCATE (bandMatches(database%noBands))
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"bandMatches")
    ALLOCATE (switchMatches(database%noSwitches))
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"switchMatches")
    ALLOCATE (spectrometerMatches(database%noSpectrometers))
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"spectrometerMatches")

    radiometerMatches=.TRUE.
    bandMatches=.TRUE.
    switchMatches=.TRUE.
    spectrometerMatches=.TRUE.

    IF (PRESENT(noCopy)) THEN
       useNoCopy=noCopy
    ELSE
       useNoCopy=.FALSE.
    ENDIF

    mostAdvancedField=0
    radiometerIndex=0
    bandIndex=0
    switchIndex=0
    spectrometerIndex=0
    channelString=""

    remainder=Capitalize(request)

    ! We go through and parse the string piece by piece
    
    ParseMLSSignalRequestWordLoop: DO
       CALL SplitWords(remainder,field,newRemainder,delimiter='.')
       IF (TRIM(field)=="") EXIT ParseMLSSignalRequestWordLoop
       keyChar=field(1:1)
       remainder=newRemainder
       SELECT CASE (keyChar)
       CASE ('R')
          IF (mostAdvancedField /= 0) &
               & CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & OutOfOrder//TRIM(request))
          CALL ParseRadiometerRequest(field,radiometerMatches)
          mostAdvancedField=1
       CASE ('B')
          IF (mostAdvancedField > 2) &
               & CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & OutOfOrder//TRIM(request))
          CALL ParseBandRequest(field,bandMatches, &
               & upperLower=upperLower)
          mostAdvancedField=2
       CASE ('S')
          IF (mostAdvancedField > 3) &
               & CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & OutOfOrder//TRIM(request))
          CALL ParseSwitchRequest(field,switchMatches)
          mostAdvancedField=3
       CASE ('C')
          channelString=field
          mostAdvancedField=5
       CASE DEFAULT             ! Must be spectrometer or wrong
          IF (mostAdvancedField > 4) &
               & CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & OutOfOrder//TRIM(request))
          CALL ParseSpectrometerRequest(field,spectrometerMatches)
          mostAdvancedField=4
       END SELECT
    END DO ParseMLSSignalRequestWordLoop

    ! Now we know what they've asked for look for complete matches

    ALLOCATE (allMatch(database%noValidSignals))
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"allMatch")
    DO signal=1,database%noValidSignals

       allMatch(signal)= &
            & radiometerMatches(database%validSignals(signal)% &
            &    radiometerIndex) .AND. &
            & bandMatches(database%validSignals(signal)%bandIndex) .AND. &
            & switchMatches(database%validSignals(signal)%switchIndex) .AND. &
            & spectrometerMatches(database%validSignals(signal)% &
            &    spectrometerIndex)
    ENDDO
       
    noMatches=COUNT(allMatch)
    IF (noMatches==0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "No matching signal:"//TRIM(request))

    ! Now create the result, if it's alread allocated that's an error
    
    IF (ASSOCIATED(signals)) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Signals already allocated")
    ALLOCATE(signals(noMatches),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"signals")

    signals=PACK(database%validSignals,allMatch)
    signals%upperLower=upperLower
    signals%notCopy=useNoCopy

    IF (.NOT. useNoCopy) CALL TurnMLSChannelInfoIntoCopy(signals)

    ! At this point we deal with the channel specifier

    IF (TRIM(channelString) /= "") THEN
       IF (useNoCopy) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Cannot use no Copy with specific channels")

       ! Also can't request specific channels when spanning more than one
       ! spectrometer type.

       IF (noMatches>1) THEN
          uniqueSpectrometerFamilyRequest=.TRUE.
          spectrometerFamily=signals(1)%spectrometerFamily
          DO signal=2,noMatches
             IF (spectrometerFamily /= signals(signal)%spectrometerFamily) &
                  & uniqueSpectrometerFamilyRequest=.FALSE.
          ENDDO
          IF (.NOT. uniqueSpectrometerFamilyRequest) CALL MLSMessage( &
               & MLSMSG_Error,ModuleName, &
               & "Cannot give specific channels accross multiple "// &
               & "spectrometer families")
       ENDIF

       CALL ParseChannelRequest(channelString, &
            & signals(1)%firstchannelInBand, signals(1)%lastChannelInBand, &
            & channelIncluded)
       DO signal=1,noMatches
          signals(signal)%channelIncluded=channelIncluded
          signals(signal)%noChannelsIncluded=COUNT(channelIncluded)
       ENDDO
       DEALLOCATE(channelIncluded, stat=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"channelIncluded")          
    ENDIF

    DEALLOCATE (radiometerMatches, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"radiometerMatches")
    DEALLOCATE (bandMatches, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"bandMatches")
    DEALLOCATE (switchMatches, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"switchMatches")
    DEALLOCATE (spectrometerMatches, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"spectrometerMatches")
    DEALLOCATE (allMatch, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"allMatch")
  END SUBROUTINE ParseMLSSignalRequest

  ! -------------------------------------------------------------------

  ! This routine destroys a signal or signal array as created by the
  ! ParseMLSSignalRequest routine.

  SUBROUTINE DestroyMLSSignalsInfo(signals,noError)

    ! Dummy arguments
    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signals
    LOGICAL, INTENT(IN), OPTIONAL :: noError

    ! Local variables
    INTEGER :: signal, status
    LOGICAL :: useNoError

    ! Executable code

    IF (PRESENT(noError)) THEN
       useNoError=noError
    ELSE
       useNoError=.FALSE.
    ENDIF

    IF (ASSOCIATED(signals)) THEN
       DO signal=1,SIZE(signals)
          IF (.NOT. signals(signal)%notCopy) THEN
             DEALLOCATE(signals(signal)%channelIncluded, stat=status)
             IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                 & MLSMSG_DeAllocate//"signals(signal)%channelIncluded")
             DEALLOCATE(signals(signal)%channelPosition, stat=status)
             IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                 & MLSMSG_DeAllocate//"signals(signal)%channelPosition")
             DEALLOCATE(signals(signal)%channelWidth, stat=status)
             IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                 & MLSMSG_DeAllocate//"signals(signal)%channelWidth")
          ENDIF
       END DO
       DEALLOCATE (signals, stat=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"signals")
    ELSE
       IF (.NOT. noError) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & "This signal not allocted")
    ENDIF
  END SUBROUTINE DestroyMLSSignalsInfo
 
  ! -------------------------------------------------------------------

  ! This subroutine combines two sets of signal lists in a set union type of
  ! operation.  Thus duplication is avoided. Note that the result must be
  ! undefined on entry. Also note, we don't consider duplication within a or b

  SUBROUTINE UnionMLSSignalsInfo(signalsA,signalsB,signalsUnion)

    ! Dummy arguments
    TYPE (MLSSignal_T), DIMENSION(:), INTENT(IN) :: signalsA
    TYPE (MLSSignal_T), DIMENSION(:), INTENT(IN) :: signalsB
    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signalsUnion

    ! Local variables
    INTEGER :: sizeA,sizeB,noSignalsInUnion
    LOGICAL, DIMENSION(:), ALLOCATABLE :: presentInA ! (sizeB)
    INTEGER :: status
    INTEGER :: noRepeats

    INTEGER :: indexA,indexB,indexUnion

    ! Executable code
    ! First work out which are the unique signals

    sizeA=SIZE(signalsA)
    sizeB=SIZE(signalsB)

    ALLOCATE (presentInA(sizeB),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"presentInA")
    presentInA=.FALSE.

    DO indexA=1,sizeA
       presentInA=presentInA .OR. (signalsA(indexA)%signalDatabaseIndex == &
            & signalsB%signalDatabaseIndex)
    ENDDO

    noRepeats=COUNT(presentInA)
    noSignalsInUnion=sizeA+sizeB-noRepeats

    ! Now allocate the result

    IF (ASSOCIATED(signalsUnion)) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Result already allocated")
    ALLOCATE (signalsUnion(noSignalsInUnion),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & MLSMSG_Allocate//"signalsUnion")

    ! Now fill up the result

    signalsUnion(1:sizeA)=signalsA
    IF (noRepeats < sizeB) signalsUnion(sizeA+1:noSignalsInUnion) = &
         & PACK(signalsB,(.NOT. presentInA))

    ! Copy channel info over rather than point to it

    CALL TurnMLSChannelInfoIntoCopy(signalsUnion)

    ! Now we need to union any channel information. We do this using a loop,
    ! but we only need to go to sizeA

    DO indexUnion=1,sizeA
       presentInA=(signalsUnion(indexUnion)%signalDatabaseIndex == &
            & signalsB%signalDatabaseIndex)
       IF (COUNT(presentInA) == 1) THEN
          indexB=1
          IndexBHunt: DO
             IF (presentInA(indexB)) EXIT IndexBHunt
             indexB=indexB+1
          END DO IndexBHunt
          signalsUnion(indexUnion)%channelIncluded= &
               & signalsUnion(indexUnion)%channelIncluded .OR. &
               & signalsB(indexB)%channelIncluded
       ENDIF
    ENDDO

    DEALLOCATE (presentInA, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"presentInA") 

  END SUBROUTINE UnionMLSSignalsInfo

  ! -------------------------------------------------------------------
  
  ! This subroutine uses the one above to concatenate two lists of signals,
  ! adding signalsB into signalsA. Note that this also allows signalsA to be
  ! unallocated, in which case signalsB is just copied.

  SUBROUTINE ConcatenateMLSSignalsInfo(signalsA,signalsB)

    ! Dummy arguments
    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signalsA
    TYPE (MLSSignal_T), DIMENSION(:), INTENT(IN) :: signalsB

    ! Local variables
    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: tempSignals=>null()
    INTEGER :: status

    ! Executable code

    IF (.NOT. ASSOCIATED(signalsA)) THEN
       ALLOCATE(signalsA(LBOUND(signalsB,1):UBOUND(signalsB,1)),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"signalsA")
       signalsA=signalsB
       CALL TurnMLSChannelInfoIntoCopy(signalsA)
    ELSE
       CALL UnionMLSSignalsInfo(signalsA,signalsB,tempSignals)
       CALL DestroyMLSSignalsInfo(signalsA)
       signalsA=>tempSignals
    ENDIF
  END SUBROUTINE ConcatenateMLSSignalsInfo

  ! -------------------------------------------------------------------

  ! This subroutine is like the union one above, but this time provides the
  ! intersection.

  SUBROUTINE IntersectionMLSSignalsInfo(signalsA,signalsB,signalsIntersection)

    ! Dummy arguments
    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signalsA
    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signalsB
    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signalsIntersection

    ! Local variables
    INTEGER :: indexA,indexB,sizeA,sizeB,status
    LOGICAL, DIMENSION(:), POINTER :: presentInA,channelIncluded

    ! Executable code

    sizeA=SIZE(signalsA)
    sizeB=SIZE(signalsB)
    ALLOCATE (presentInA(sizeB),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & MLSMSG_Allocate//"presentInA")

    IF (ASSOCIATED(signalsIntersection)) CALL MLSMessage( &
         & MLSMSG_Error,ModuleName,"Result already allocated")

    ! Unlike the union routine this one has to be a little more do loop
    ! orientated, due to the use of channels

    DO indexA=1,sizeA
       presentInA= (signalsA(indexA)%signalDatabaseIndex == &
            & signalsB%signalDatabaseIndex)
       IF (COUNT(presentInA) /= 0) THEN
          indexB=1
          IndexBHunt: DO
             IF (presentInA(indexB)) EXIT IndexBHunt
             indexB=indexB+1
          END DO IndexBHunt

          ! Now look at the channels

          ALLOCATE(channelIncluded( &
               & signalsA(indexA)%firstChannelInBand: &
               & signalsA(indexA)%lastChannelInBand),STAT=status)
          IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & "Allocation failed for channelIncluded")
          
          channelIncluded=signalsA(indexA)%channelIncluded .AND. &
               & signalsB(indexB)%channelIncluded

          IF (COUNT(channelIncluded) /= 0) THEN
             CALL ConcatenateMLSSignalsInfo(signalsIntersection, &
                  & signalsA(indexA:indexA))
             signalsIntersection(UBOUND(signalsIntersection,1))% &
                  & channelIncluded = channelIncluded
          ENDIF
          DEALLOCATE (channelIncluded, stat=status)       
          IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
              & MLSMSG_DeAllocate//"channelIncluded")
       ENDIF
    ENDDO

    DEALLOCATE (presentInA, stat=status)       
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"presentInA")
  END SUBROUTINE IntersectionMLSSignalsInfo

  ! -------------------------------------------------------------------

  ! This subroutine gives the full name for an MLS signal

  SUBROUTINE GetFullMLSSignalName(signal,fullName,rbOnly,showChannel)
    
    ! Dummy arguments
    TYPE (MLSSignal_T), INTENT(IN) :: signal
    CHARACTER (LEN=*), INTENT(OUT) :: fullName
    LOGICAL, OPTIONAL, INTENT(IN) :: rbOnly ! Only show radiometer and band
    LOGICAL, OPTIONAL, INTENT(IN) :: showChannel ! Only show radiometer and band

    ! Local variables
    LOGICAL :: useRBOnly
    LOGICAL :: useShowChannel
    LOGICAL :: previousChanIncluded,thisChanIncluded
    INTEGER :: channel,firstChannelInRange,rangeLen
    LOGICAL :: firstEntry
    CHARACTER (LEN=32) :: word

    ! Executable code

    fullName=TRIM(signal%radiometerName)//"."//TRIM(signal%bandName)

    IF (PRESENT(rbOnly)) THEN
       useRBOnly=rbOnly
    ELSE
       useRBOnly=.FALSE.
    ENDIF

    IF (PRESENT(showChannel)) THEN
       useShowChannel=showChannel
    ELSE
       useShowChannel=.FALSE.
    ENDIF

    IF (.NOT. useRBOnly) THEN
       fullName=TRIM(fullName)//"."//TRIM(signal%switch)//"."// &
            TRIM(signal%spectrometerName)
    ENDIF

    IF (useShowChannel) THEN
       fullName=TRIM(fullName)//".C["
       firstEntry=.TRUE.
       channel=signal%firstChannelInBand-1
       previousChanIncluded=.FALSE.
       ShowChannelLoop: DO
          channel=channel+1
          IF (channel<=signal%noChannelsInBand) THEN
             thisChanIncluded=signal%channelIncluded(channel)
          ELSE
             thisChanIncluded=.FALSE.
          ENDIF

          IF (thisChanIncluded.NEQV.previousChanIncluded) THEN
             ! We have a transition of some kind
             IF (thisChanIncluded) THEN ! Begining of a new range
                firstChannelInRange=channel
                WRITE (UNIT=word,FMT=*) firstChannelInRange
                IF (.NOT. firstEntry) fullName=TRIM(fullName)//"+"
                firstEntry=.FALSE.
                fullName=TRIM(fullName)//TRIM(ADJUSTL(word))
             ELSE               ! The end of a range or of array
                rangeLen=channel-firstChannelInRange
                IF (rangeLen>1) THEN
                   WRITE (UNIT=word,FMT=*) channel-1
                   IF (rangeLen==2) THEN
                      fullName=TRIM(fullName)//"+"//TRIM(ADJUSTL(word))
                   ELSE
                      fullName=TRIM(fullName)//":"//TRIM(ADJUSTL(word))
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
          IF (channel==signal%lastChannelInBand+1) EXIT ShowChannelLoop
          previousChanIncluded=thisChanIncluded
       END DO ShowChannelLoop
       fullName=TRIM(fullName)//"]"
    ENDIF

    ! We might put channel stuff here later

  END SUBROUTINE GetFullMLSSignalName

  ! -------------------------------------------------------------------

  ! This routine reads the signals database file and fills the data structure
  ! The routine doesn't do very much error checking as this file is assumed not
  ! to change very often.

  SUBROUTINE ReadSignalsDatabase(unit)

    ! Arguments and result

    INTEGER, INTENT(IN) :: unit

    ! Local parameters ------------------------------

    INTEGER, PARAMETER :: SDBLineLen=132
    CHARACTER (LEN=*), PARAMETER :: EOFMessage= &
         & "Unexpected EOF on signals database file"

    ! Some low level variables ----------------------

    CHARACTER (LEN=SDBLineLen) :: line, first, last, rest
    LOGICAL :: eof
    INTEGER :: no               ! Temporary array index
    INTEGER :: signal           ! Loop counters
    INTEGER :: radiometer, band, switch, spectrometer, channel ! Loop counters
    INTEGER :: firstChannel, lastChannel, noChannels
    INTEGER :: spectrometerFamily ! Another loop counter
    INTEGER :: wordLen          ! Length of word
    INTEGER :: hasModifier      ! Flag
    INTEGER :: status           ! From allocate
    INTEGER :: index            ! General array index

    INTEGER :: evenNo           ! For `even' channels
    REAL(r8) :: evenStart,evenSpacing,evenWidth ! For `even' channels

    TYPE (MLSSignal_T), DIMENSION(:), POINTER :: tempSignal=>null()

    ! These variables are intermediate arrays to allow our database to `grow'

    TYPE (SDBSpectrometerFamilyInfo_T), DIMENSION(:), POINTER :: &
         & tempSpectrometerFamilyInfo
    CHARACTER (LEN=SDBLineLen), DIMENSION(:), POINTER :: &
         validSignalNames, tempValidSignalNames

    ! These strings are a breakdown of the valid signals strings

    CHARACTER (LEN=NameLen), DIMENSION(:), ALLOCATABLE :: validRadiometer, &
         & validBand, validSwitch, validSpectrometer
    CHARACTER (LEN=NameLen), DIMENSION(:), POINTER :: radiometerNames, &
         & bandNames, switchNames, spectrometerNames

    CHARACTER (LEN=1), DIMENSION(:), ALLOCATABLE :: &
         & spectrometerFamilyChars


    ! Executable code ----------------------------------------------------

    ! The first section in the signals file describes the various types
    ! of spectrometer there can be.  First, we read a line and expect it to be
    ! `spectrometers'
    INTEGER, PARAMETER :: MLSInstrumentNoModules=2
    CHARACTER (LEN=3),  DIMENSION(MLSInstrumentNoModules) :: &
       & MLSInstrumentModuleNames= (/ &
       & "GHz", &
       & "THz"/)


    MLSInstrumentModuleNames(1)='GHz'
    MLSInstrumentModuleNames(2)='THz'
    !for God's sake!

    CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
    IF (eof.OR.(Capitalize(TRIM(line))/="SPECTROMETERS")) &
         & CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "SPECTROMETERS expected in signals database file")

    ! Now we go through and read the family descriptions

    database%noSpectrometerFamilies=0
    SpectrometerFamilyLoop: DO
       ! Read the line
       CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
       IF (eof) CALL MLSMessage(MLSMSG_Error,ModuleName,EOFMessage)
       IF (Capitalize(TRIM(line))=="END") EXIT SpectrometerFamilyLoop

       database%noSpectrometerFamilies=database%noSpectrometerFamilies+1
       no=database%noSpectrometerFamilies
       IF (no==1) THEN
          ALLOCATE(database%spectrometerFamilyInfo(no),STAT=status)
          IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & MLSMSG_Allocate//"database%spectrometerFamilyInfo")
       ELSE
          ALLOCATE(tempSpectrometerFamilyInfo(no),STAT=status)
          IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & MLSMSG_Allocate//"tempSpectrometerFamilyInfo")
          tempSpectrometerFamilyInfo(1:no-1)=database%spectrometerFamilyInfo
          DEALLOCATE(database%spectrometerFamilyInfo, stat=status)
          IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             & MLSMSG_DeAllocate//"database%spectrometerFamilyInfo")

          database%spectrometerFamilyInfo=>tempSpectrometerFamilyInfo
       ENDIF

       ! Now parse the line

       CALL SplitWords(line,first,rest,last,threeWay=.TRUE.,delimiter=" ")
       database%spectrometerFamilyInfo(no)%name=TRIM(first)
       READ (UNIT=rest,FMT=*) database%spectrometerFamilyInfo(no)%firstChannel
       READ (UNIT=last,FMT=*) database%spectrometerFamilyInfo(no)%lastChannel
       database%spectrometerFamilyInfo(no)%noChannels= &
            & database%spectrometerFamilyInfo(no)%lastChannel- &
            & database%spectrometerFamilyInfo(no)%firstChannel+1
    END DO SpectrometerFamilyLoop

    ! The next section in the file is the list of all the valid signals
    ! We'll just read these in for the moment, and parse them in the next
    ! section

    CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
    IF (eof) CALL MLSMessage(MLSMSG_Error,ModuleName,EOFMessage)
    IF (Capitalize(TRIM(line))/="VALID SIGNALS") CALL MLSMessage( &
         & MLSMSG_Error,ModuleName,"Valid signals expected")

    database%noValidSignals=0
    ValidSignalsReadingLoop: DO
       ! Read the line
       CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
       IF (eof) CALL MLSMessage(MLSMSG_Error,ModuleName,EOFMessage)
       IF (Capitalize(TRIM(line))=="END") EXIT ValidSignalsReadingLoop

       database%noValidSignals=database%noValidSignals+1
       no=database%noValidSignals
       IF (no==1) THEN
          ALLOCATE(validSignalNames(no),STAT=status)
          IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & MLSMSG_Allocate//"validSignalNames")
       ELSE
          ALLOCATE(tempValidSignalNames(no),STAT=status)
          IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & MLSMSG_Allocate//"tempValidSignalNames")
          tempValidSignalNames(1:no-1)=validSignalNames
          DEALLOCATE(validSignalNames, stat=status)
          IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             & MLSMSG_DeAllocate//"validSignalNames")
          validSignalNames=>tempValidSignalNames
       ENDIF
       validSignalNames(no)=line
    END DO ValidSignalsReadingLoop

    ! The remainder of the file talks about lo frequencies etc.  We'll handle
    ! all that stuff in a minute.  First we'll take apart that information we
    ! got

    ALLOCATE(validRadiometer(database%noValidSignals),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         & "validRadiometer")
    ALLOCATE(validBand(database%noValidSignals),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         & "validBand")
    ALLOCATE(validSwitch(database%noValidSignals),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         & "validSwitch")
    ALLOCATE(validSpectrometer(database%noValidSignals),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         & "validSpectrometer")

    ALLOCATE(radiometerNames(database%noValidSignals),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         & "radiometerNames")
    ALLOCATE(bandNames(database%noValidSignals),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         & "bandNames")
    ALLOCATE(switchNames(database%noValidSignals),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         & "switchNames")
    ALLOCATE(spectrometerNames(database%noValidSignals),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         & "spectrometerNames")

    DO signal=1,database%noValidSignals
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

    ALLOCATE(database%radiometerInfo(database%noRadiometers),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate // "database%radiometerInfo")

    DO radiometer=1,database%noRadiometers
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

    ALLOCATE(database%bandInfo(database%noBands),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"database%bandInfo")

    DO band=1,database%noBands
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

    ALLOCATE ( database%switches(database%noSwitches), STAT=status )
    IF ( status /= 0 ) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"database%switches")
    database%switches=switchNames(1:database%noSwitches)

    ! Finally we do the spectrometers

    ALLOCATE (database%spectrometerInfo(database%noSpectrometers), &
         &  STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"database%spectrometerInfo")
    DO spectrometer=1,database%noSpectrometers
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

    ALLOCATE (spectrometerFamilyChars(database%noSpectrometerFamilies))

    DO spectrometerFamily=1,database%noSpectrometerFamilies
       spectrometerFamilyChars(spectrometerFamily)= &
            & database%spectrometerFamilyInfo(spectrometerFamily)%name(1:1)
    END DO

    DO spectrometer=1,database%noSpectrometers
       database%spectrometerInfo(spectrometer)%familyIndex= &
            & LinearSearchStringArray(database%spectrometerFamilyInfo%name, &
            & database%spectrometerInfo(spectrometer)%fullFamily)
    END DO

    DO band=1,database%noBands
       database%bandInfo(band)%spectrometerFamilyIndex= &
            & LinearSearchStringArray(spectrometerFamilyChars, &
            & database%bandInfo(band)%spectrometerFamily)
    END DO

    ! That's basically the valid signals database filled up, no we go on and
    ! read the rest of the file which contains frequency information etc.

    ! The next section descusses the radiometer frequencies

    CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
    IF (eof.OR.(Capitalize(TRIM(line))/="RADIOMETER FREQUENCIES")) &
         & CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "RADIOMETER FREQUENCIES expected in database")

    ! Now this section is just a list of names followed by lo's

    database%radiometerInfo%lo=0.0D0

    RadiometerFrequencyLoop: DO
       CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
       IF (eof) CALL MLSMessage(MLSMSG_Error,ModuleName,EOFMessage)
       IF (Capitalize(TRIM(line))=="END") EXIT RadiometerFrequencyLoop
       CALL SplitWords(line,first,rest,last,delimiter=' ',threeWay=.TRUE.)
       radiometer=LinearSearchStringArray(database%radiometerInfo%name, &
            & first, caseInsensitive=.TRUE.)
       IF (radiometer == 0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
            & "No such radiometer: "//TRIM(first))
       READ (unit=rest,FMT=*) database%radiometerInfo(radiometer)%lo
       database%radiometerInfo(radiometer)%instrumentModule=&
!           & LinearSearchStringArray(MLSInstrumentModuleNames(:),last,&
            & LinearSearchStringArray(MLSInstrumentModuleNames,last,&
            & caseInsensitive=.TRUE.)

       IF (database%radiometerInfo(radiometer)%instrumentModule==&
            & L_None) CALL MLSMessage(MLSMSG_Error,&
            & ModuleName,"Unrecognised instrument module: "//last)
    END DO RadiometerFrequencyLoop

    IF (MINVAL(database%radiometerInfo%lo) <= 0.0D0) CALL MLSMessage(&
         & MLSMSG_Error,ModuleName,"Not all radiometer LOs assigned")

    ! Now the next section is a set of similar information for the bands

    CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
    IF (eof.OR.(Capitalize(TRIM(line))/="BAND FREQUENCIES")) &
         & CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "BAND FREQUENCIES expected in database")

    ! Now this section is just a list of names followed by center frequencies

    database%bandInfo%centerFreqIF=0.0D0

    BandFrequencyLoop: DO
       CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
       IF (eof) CALL MLSMessage(MLSMSG_Error,ModuleName,EOFMessage)
       IF (Capitalize(TRIM(line))=="END") EXIT BandFrequencyLoop
       CALL SplitWords(line,first,rest,delimiter=" ")
       band=LinearSearchStringArray(database%bandInfo%name, &
            & first, caseInsensitive=.TRUE.)
       IF (band == 0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
            & "No such band: "//first)
       READ (unit=rest,FMT=*) database%bandInfo(band)%centerFreqIF
    END DO BandFrequencyLoop

    ! Don't check here for all allocated as the WF type filters are listed
    ! seperately.

    ! Now we have the section giving the frequencies and widths of all the
    ! spectrometer channels.

    CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
    IF (eof.OR.(Capitalize(TRIM(line))/="SPECTROMETER FREQUENCIES")) &
         & CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "SPECTROMETER FREQUENCIES expected in database")

    SpectrometerFrequencyLoop: DO
       CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
       IF (eof) CALL MLSMessage(MLSMSG_Error,ModuleName,EOFMessage)
       IF (Capitalize(TRIM(line))=="END") EXIT SpectrometerFrequencyLoop

       ! Each spectrometer family has a line describing it with a set of
       ! instructions as to how to assign the channels

       CALL SplitWords(line,first,rest,delimiter=" ")
       index=LinearSearchStringArray( &
            & database%spectrometerFamilyInfo%name,first, &
            & caseInsensitive=.TRUE.)
       IF (index == 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & "No such spectrometer family: "//first)

       ! Here we'll allocate the arrays for size and posiiton.  Note for the WF
       ! filters this is unneccessary, and we'll deallocate them later.
       ! But to make the code more readable I'm allocating them first.

       ALLOCATE (database%spectrometerFamilyInfo(index)%position(&
            & database%spectrometerFamilyInfo(index)%firstChannel:&
            & database%spectrometerFamilyInfo(index)%lastChannel),&
            & STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//first//" position")

       ALLOCATE (database%spectrometerFamilyInfo(index)%width(&
            & database%spectrometerFamilyInfo(index)%firstChannel:&
            & database%spectrometerFamilyInfo(index)%lastChannel),&
            & STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//first//" width")

       database%spectrometerFamilyInfo(index)%individual=.FALSE.

       SELECT CASE (Capitalize(TRIM(rest)))
       CASE ('LIST') ! Positions and widths given on next lines
          READ (UNIT=unit, FMT=*) &
               & database%spectrometerFamilyInfo(index)%position
          READ (UNIT=unit, FMT=*) &
               & database%spectrometerFamilyInfo(index)%width

       CASE ('EVEN') ! No chans, start Freq, spacing and widths
          READ (UNIT=unit, FMT=*) evenNo,evenStart,evenSpacing,evenWidth
          IF (evenNo /= database%spectrometerFamilyInfo(index)%noChannels) &
               & CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & "Wrong number of channels for "//first)
          ! Loop and fill information
          database%spectrometerFamilyInfo(index)%width=evenWidth
          DO evenNo=database%spectrometerFamilyInfo(index)%firstChannel, &
               & database%spectrometerFamilyInfo(index)%lastChannel
             database%spectrometerFamilyInfo(index)%position(evenNo)= &
                  & evenStart+evenNo*evenSpacing
          END DO

       CASE ('INDIVIDUAL')
          database%spectrometerFamilyInfo(index)%individual=.TRUE.
          DEALLOCATE(database%spectrometerFamilyInfo(index)%position, stat=status)
          IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"database%spectrometerFamilyInfo(index)%position") 
          DEALLOCATE(database%spectrometerFamilyInfo(index)%width, stat=status)
          IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"database%spectrometerFamilyInfo(index)%width")

       CASE DEFAULT 
          CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & "Unrecognized spectrometer family: "//first)
       END SELECT
    END DO SpectrometerFrequencyLoop

    ! The final section of the file is the section dealing with specific
    ! channels. Here, having built up most of the database, we can find the
    ! relevant data using our standard parse routines.

    ! Now we have read all the gory information from our database
    ! we need to fill up our valid signals database

    ALLOCATE (database%validSignals(database%noValidSignals),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & MLSMSG_Allocate//"database%validSignals")

    DO signal=1,database%noValidSignals

       database%validSignals(signal)%signalDatabaseIndex=signal

       ! Deal with radiometer info first.

       radiometer=LinearSearchStringArray( &
            & database%radiometerInfo%name,validRadiometer(signal))
       database%validSignals(signal)%radiometerIndex=radiometer
       database%validSignals(signal)%instrumentModule=&
            & database%radiometerInfo(radiometer)%instrumentModule
       database%validSignals(signal)%radiometerName= &
            & database%radiometerInfo(radiometer)%name
       database%validSignals(signal)%radiometerPrefix= &
            & database%radiometerInfo(radiometer)%prefix
       database%validSignals(signal)%radiometerSuffix= &
            & database%radiometerInfo(radiometer)%suffix
       database%validSignals(signal)%radiometerNumber= &
            & database%radiometerInfo(radiometer)%number
       database%validSignals(signal)%radiometerModifier= &
            & database%radiometerInfo(radiometer)%modifier
       database%validSignals(signal)%lo= &
            & database%radiometerInfo(radiometer)%lo      

       ! Now the band info

       band=LinearSearchStringArray(database%bandInfo%name,validBand(signal))
       database%validSignals(signal)%bandIndex=band
       database%validSignals(signal)%bandName=database%bandInfo(band)%name
       database%validSignals(signal)%bandSuffix=database%bandInfo(band)%suffix
       database%validSignals(signal)%spectrometerFamily= &
            & database%bandInfo(band)%spectrometerFamily
       database%validSignals(signal)%bandCenterFreqIF= &
            & database%bandInfo(band)%centerFreqIF

       ! Now the switch info

       switch=LinearSearchStringArray(database%switches,validSwitch(signal))
       database%validSignals(signal)%switchIndex=switch
       database%validSignals(signal)%switch=database%switches(switch)

       ! Now the spectrometer info

       spectrometer=LinearSearchStringArray(database%spectrometerInfo%name, &
            & validSpectrometer(signal))
       spectrometerFamily=database%spectrometerInfo(spectrometer)%familyIndex

       database%validSignals(signal)%spectrometerIndex=spectrometer
       database%validSignals(signal)%spectrometerFamilyIndex=spectrometerFamily

       database%validSignals(signal)%spectrometerName= &
            & database%spectrometerInfo(spectrometer)%name
       database%validSignals(signal)%fullSpectrometerFamily= &
            & database%spectrometerInfo(spectrometer)%fullFamily
       database%validSignals(signal)%spectrometerNumber= &
            & database%spectrometerInfo(spectrometer)%number

       ! Now the information on the channels

       firstChannel=database% &
            & spectrometerFamilyInfo(spectrometerFamily)%firstChannel
       lastChannel=database% &
            & spectrometerFamilyInfo(spectrometerFamily)%lastChannel
       noChannels=database% &
            & spectrometerFamilyInfo(spectrometerFamily)%noChannels

       ! Now we compute the channel positions in IF space
       
       database%validSignals(signal)%firstChannelInBand=firstChannel
       database%validSignals(signal)%lastChannelInBand=lastChannel
       database%validSignals(signal)%noChannelsInBand=noChannels
       
       ALLOCATE (database%validSignals(signal)%channelPosition &
            & (firstChannel:lastChannel),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"channelPosition")

       ALLOCATE (database%validSignals(signal)%channelWidth &
            & (firstChannel:lastChannel),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"channelWidth")

       ALLOCATE (database%validSignals(signal)%channelIncluded &
            & (firstChannel:lastChannel),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"channelIncluded")

       database%validSignals(signal)%channelIncluded=.TRUE.

       DO channel=firstChannel,lastChannel
          ! We'll have to defer dealing with the individual channels
          IF (.NOT. database%spectrometerFamilyInfo(spectrometerFamily)% &
               & individual) THEN
             database%validSignals(signal)%channelPosition(channel)= &
                  & database%validSignals(signal)%bandCenterFreqIF + &
                  & database%spectrometerFamilyInfo(spectrometerFamily)% &
                  & position(channel)
             database%validSignals(signal)%channelWidth(channel)= &
                  & database%spectrometerFamilyInfo(spectrometerFamily)% &
                  & width(channel)
          ENDIF
       END DO                   ! Channel loop
    END DO                      ! Signal loop

    ! The next section of the file deals with the WF type channels.  This is
    ! actually very easy.  We simply read the spec and parse it as we've done
    ! others and store it in our database, making use of the
    ! ParseMLSSignalRequest routine which now works because we've fill up the
    ! daabase.

    CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
    IF (eof.OR.(Capitalize(TRIM(line))/="SPECIFIC CHANNELS")) &
         & CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "SPECIFIC CHANNELS expected in database")

    SpecificChannelLoop: DO
       CALL ReadCompleteLineWithoutComments(unit,line,eof=eof)
       IF (eof) CALL MLSMessage(MLSMSG_Error,ModuleName,EOFMessage)

       ! This line will either be END or a spec followed by the word 'List'
       line=Capitalize(line)
       IF (TRIM(line)=="END") EXIT SpecificChannelLoop

       CALL SplitWords(line,first,last,delimiter=" ")
       IF (TRIM(last)/="LIST") CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & "List expected for specific channel information")

       ! This will be followed by a list of positions and widths, we'll read
       ! this directory into our database by using the noCopy option.

       CALL ParseMLSSignalRequest(first,tempSignal,noCopy=.TRUE.)
       IF (SIZE(tempSignal)/=1) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & "Ambiguous spec.:"//first)

       READ (UNIT=unit,FMT=*) tempSignal(1)%channelPosition, &
            & tempSignal(1)%channelWidth

       ! Now destroy this temporary information

       CALL DestroyMLSSignalsInfo(tempSignal)
    END DO SpecificChannelLoop
       

          
    ! Now we tidy up our arrays and exit

    DEALLOCATE(validSignalNames, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"validSignalNames")
    DEALLOCATE(radiometerNames, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"radiometerNames")
    DEALLOCATE(bandNames, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"bandNames")
    DEALLOCATE(switchNames, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"switchNames")
    DEALLOCATE(spectrometerNames, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"spectrometerNames")
    DEALLOCATE(spectrometerFamilyChars, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"spectrometerFamilyChars")
  END SUBROUTINE ReadSignalsDatabase

  ! -------------------------------------------------------------------

  ! This routine deallocates all the information dealing with the signals
  ! database

  SUBROUTINE DestroySignalsDatabase

    ! Local variable

    INTEGER :: i, status

    ! Executable code

    database%noRadiometers=0
    database%noBands=0
    database%noSwitches=0
    database%noSpectrometers=0
    database%noSpectrometerFamilies=0
    database%noValidSignals=0

    DEALLOCATE (database%radiometerInfo, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"database%radiometerInfo")
    DEALLOCATE (database%bandInfo, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"database%bandInfo")
    DEALLOCATE (database%switches, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"database%switches")
    DEALLOCATE (database%spectrometerInfo, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"database%spectrometerInfo")

    DO i=1,database%noSpectrometerFamilies
       DEALLOCATE (database%spectrometerFamilyInfo(i)%position, stat=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"database%spectrometerFamilyInfo(i)%position")
       DEALLOCATE (database%spectrometerFamilyInfo(i)%width, stat=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"database%spectrometerFamilyInfo(i)%width")
    END DO
    DEALLOCATE (database%spectrometerFamilyInfo, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"database%spectrometerFamilyInfo")

    DO i=1,database%noValidSignals
       DEALLOCATE (database%validSignals(i)%channelIncluded, stat=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"database%validSignals(i)%channelIncluded")
       DEALLOCATE (database%validSignals(i)%channelPosition, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"database%validSignals(i)%channelPosition")
       DEALLOCATE (database%validSignals(i)%channelWidth, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"database%validSignals(i)%channelWidth")
    END DO
    DEALLOCATE (database%validSignals, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             &MLSMSG_DeAllocate//"database%validSignals")

  END SUBROUTINE DestroySignalsDatabase

  ! -------------------------------------------------------------------

  ! This routine returns an array of the MLS radiometer names from the database

  SUBROUTINE GetMLSRadiometerNames(names)

    ! Dummy arguments
    CHARACTER (LEN=NameLen), DIMENSION(:), POINTER :: names

    ! Local variables
    INTEGER :: status

    ! Executable code

    ALLOCATE(names(database%noRadiometers),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         & "names")

    names=database%radiometerInfo%name
  END SUBROUTINE GetMLSRadiometerNames

  ! -------------------------------------------------------------------

  ! This routine returns an array of names of bands in MLS from the database.

  SUBROUTINE GetMLSBandNames(names)

    ! Dummy arguments
    CHARACTER (LEN=NameLen), DIMENSION(:), POINTER :: names

    ! Local variables
    INTEGER :: status

    ! Executable code

    ALLOCATE(names(database%noBands),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         & "names")

    names=database%bandInfo%name
  END SUBROUTINE GetMLSBandNames

!=============================================================================
END MODULE MLSSignalNomenclature
!=============================================================================

!
! $Log$
! Revision 2.3  2001/03/30 19:51:49  pwagner
! *** empty log message ***
!
! Revision 2.1  2001/02/09 00:38:56  livesey
! Various changes
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.21  2000/06/30 00:59:38  lungu
! Initialized tempsignal, tempsignals=>null()
! Changed switchNames to switchNames(1:database%noSwitches)
! so that a "shape error" is avoided.
!
