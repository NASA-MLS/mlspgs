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
module MLSSignalNomenclature    ! Dealing with MLS rad.band etc. specifiers
!=============================================================================

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_test, &
    & Test_Allocate, Test_Deallocate
  use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
  use MLSCommon, only: NameLen
  use MLSKinds, only: R8
  use MLSStrings, only: Capitalize, LinearSearchStringArray, &
    & ReadCompleteLineWithoutComments, splitWords
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error

  implicit none
  private

  public :: MLSSignal_T, ParseMLSSignalRequest, DestroyMLSSignalsInfo,&
       & ConcatenateMLSSignalsInfo, UnionMLSSignalsInfo, &
       & IntersectionMLSSignalsInfo, GetFullMLSSignalName, &
       & ReadSignalsDatabase, DestroySignalsDatabase, &
       & GetMLSRadiometerNames, GetMLSBandNames

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

  ! The same as 2.1 before being removed (see log below)

  ! This module contains all the functionality required for dealing with the
  ! standard MLS signal nomenclature.  The instrument configuration is read
  ! from a file (e.g. emls-signals.dat). The user can request information about
  ! a signal or set of signals by name.

  ! This datatype describes a valid radiometer/switch/band/spectrometer
  ! combination.  It is the one most used by calling code.  For this reason I
  ! have put it first.  Some of the information given is duplicated in later
  ! structures but they are typically only internally used by this module.

  type MLSSignal_T

     ! First we have a set of indices into arrays.  These are used a lot in the
     ! module.  But the calling code may also find them useful for issues such
     ! as: Does signal A come from the same radiometer as signal B?

     integer :: signalDatabaseIndex ! Index into rad.band.s.spec combination
     integer :: instrumentModule ! Module in the instrument (EMLS:1=GHz, 2=THz)
     integer :: radiometerIndex ! Index into array of radiometerInfos
     integer :: bandIndex       ! Index into array of bandInfos
     integer :: switchIndex     ! Index into array of switch names
     integer :: spectrometerIndex ! Index into spectrometer array
     integer :: upperLower      ! -1=Lower, 1=Upper, 0=Folded
     logical :: notCopy         ! If set POINTER arrays below are not copies

     ! For an explanation of the notCopy parameter see the
     ! ParseMLSSignalRequest routine.

     ! Now we have detailed information on the radiometer

     character (len=NameLen) :: radiometerName   ! e.g. R1A:118
     character (len=NameLen) :: radiometerPrefix ! e.g. R1A
     character (len=NameLen) :: radiometerSuffix ! e.g. 118
     integer :: radiometerNumber                    ! e.g. 1
     character (len=1) :: radiometerModifier        ! e.g. A
     real(r8) :: lo                         ! MHz

     ! Now we have detailed information on the band

     character (len=NameLen) :: bandName    ! e.g. B1F:PT
     character (len=NameLen) :: bandSuffix  ! e.g. PT
     character (len=1) :: spectrometerFamily   ! e.g. F
     real(r8) :: bandCenterFreqIF      ! MHz

     ! Now we have information on the switch

     character (len=NameLen) :: switch

     ! Now information on the spectrometers

     character (len=NameLen) :: spectrometerName
     character (len=NameLen) :: fullSpectrometerFamily ! e.g. FB25
     integer :: spectrometerFamilyIndex ! Index into array of spec. fams.
     integer :: spectrometerNumber ! Note count from one for this.

     ! Now information the channels in the whole band

     integer :: firstChannelInBand
     integer :: lastChannelInBand
     integer :: noChannelsInBand

     ! Now this particular signal may be a subset of all the channels so
     ! detail what we have.

     integer :: noChannelsIncluded
     logical, dimension(:), pointer :: channelIncluded
     real(r8), dimension(:), pointer :: channelPosition ! i.f. space
     real(r8), dimension(:), pointer :: channelWidth ! i.f. space

  end type MLSSignal_T

  ! ----------------------------------------------------------------------

  ! The remaining datatypes are somewhat more private.

  ! This datatype describes a radiometer
  type SDBRadiometerInfo_T
     character (len=nameLen) :: name   ! e.g. R1A:118
     character (len=nameLen) :: prefix ! e.g. R1A
     character (len=nameLen) :: suffix ! e.g. 118
     integer :: number                     ! e.g. 1
     character(len=1) :: modifier          ! e.g. A/B or H/V for R5 (emls)
     real(r8) :: lo                ! Local oscillator /MHz
     integer :: instrumentModule ! Module in instrument GHz/THz
  end type SDBRadiometerInfo_T

  ! This datatype describes a band
  type SDBBandInfo_T
     character (len=NameLen) :: name      ! e.g. B1F:PT
     character (len=NameLen) :: suffix    ! e.g. PT
     integer :: number                       ! e.g. 1
     character (len=1) :: spectrometerFamily ! e.g. F
     integer :: spectrometerFamilyIndex      ! Index into array of next type
     real(r8) :: centerFreqIF        ! Center i.f. frequency (MHz)
  end type SDBBandInfo_T

  ! This datatype describes a spectrometer family
  type SDBSpectrometerFamilyInfo_T
     character (len=NameLen) :: name ! Name of family e.g. FB25
     integer :: noSpectrometers ! Number of spectrometers in this family
     integer :: noChannels      ! Number of channels in family
     integer :: firstChannel    ! First channel number (e.g. 1)
     integer :: lastChannel     ! Last channel number (e.g. 25)
     logical :: individual      ! If set have discrete freqs. e.g. wf4 series

     ! These two arrays are the position and width of the channels wrt. the if
     ! The arrays are actually dimensioned firstChannel:lastChannel.
     ! Units are MHz

     real(r8), dimension(:), pointer :: position
     real(r8), dimension(:), pointer :: width
  end type SDBSpectrometerFamilyInfo_T

  ! This small datatype describes a spectrometer
  type SDBSpectrometerInfo_T
     character (len=NameLen) :: name ! Name of spectrometer
     character (len=NameLen) :: fullFamily ! Full name of family eg FB25
     character (len=1) :: family ! Single character family id e.g. F
     integer :: familyIndex     ! Index into familyInfo database
     integer :: number          ! Number within family
  end type SDBSpectrometerInfo_T

  ! This datatype is an amalgam of the above and is the database that is filled
  type MLSSignalsDatabase_T
     integer :: noRadiometers   ! Including redundant etc.
     integer :: noBands         ! Accross the whole instrument
     integer :: noSwitches      ! No. vald S0, S1 etc. fields
     integer :: noSpectrometers ! No. spectrometers in whole instrument
     integer :: noSpectrometerFamilies
     integer :: noValidSignals    ! Number of valid Signal combinations

     type (SDBRadiometerInfo_T), dimension(:), pointer :: radiometerInfo => null()
          ! Actually dimensioned (noRadiometers)
     type (SDBBandInfo_T), dimension(:), pointer :: bandInfo => null()
          ! Actually dimensioned (noBands)
     character (len=NameLen), dimension(:), pointer :: switches => null()
          ! Actually dimensioned (noSwitches)
     type (SDBSpectrometerInfo_T), dimension(:), pointer :: spectrometerInfo => null()
     type (SDBSpectrometerFamilyInfo_T), dimension(:), pointer :: &
          & spectrometerFamilyInfo => null()
          ! Actually dimensioned (noSpectrometerFamilies)
     type (MLSSignal_T), dimension(:), pointer :: validSignals => null()
          ! Valid Signal combinations
          ! Actually dimensioned (noValidSignals)
  end type MLSSignalsDatabase_T

  ! Local, private variables

  type (MLSSignalsDatabase_T) :: database

contains

  ! ====================================================================

  ! -------------------------------------  ParseRadiometerRequest  -----

  subroutine ParseRadiometerRequest ( request, matches )

    ! Parse a request for a particular radiometer
    ! and return an array of flags indicating which radiometers match

    ! Dummy arguments
    character (len=*), intent(in) :: request

    logical, dimension(:) :: matches ! Result (database%noRadiometers)

    ! Local variables
    character (len=len(request)) :: prefix, suffix
    character (len=1) :: modifierRequested
    integer :: numberRequested, hasModifier, prefixLen, radiometer

    ! Executable code

    call SplitWords(request,prefix,suffix,delimiter=':')
    prefix=Capitalize(prefix)
    suffix=Capitalize(suffix)

    if (prefix(1:1)/="R") call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "R expected in radiometer specifier" )

    ! Now parse the rest of the prefix.
    ! The format of the field is R<number><modifier>, where modifier may be
    ! omitted or set to *. Also R* is valid, or just R!

    prefixLen = LEN_TRIM(prefix)
    numberRequested = 0
    modifierRequested = "*"
    hasModifier = 1

    ! First we'll look for the modifier, look at the last characater.  If it's
    ! numerical set the modifier to "*"

    if ( prefixLen>1 ) then
       modifierRequested=prefix(prefixLen:prefixLen)
       if ( (LGE(modifierRequested,"0")).AND.(LLE(modifierRequested,"9")) ) then
          hasModifier=0
          modifierRequested="*"
       end if

       ! Now we look at the numeric field if there is one.

       if ( prefixLen-hasModifier > 1 ) &
            & read (UNIT=prefix(2:prefixLen-hasModifier),FMT=*) numberRequested
    end if

    ! Now we have the requsted number (or 0 if dont care) and requested
    ! modifier (or * if dont care)

    if ( size(matches) /= database%noRadiometers) call MLSMessage ( &
         MLSMSG_Error, ModuleName, "Result is wrong size" )

    matches= &
         & ( (numberRequested == database%radiometerInfo%number) .OR. &
         &   (numberRequested == 0) ) .AND. &
         & ( (modifierRequested == database%radiometerInfo%modifier) .OR. &
         &   (modifierRequested == "*") )

    ! To match the suffix we have to be a little more old fasioned because we
    ! want to use trim.

    if ( suffix /= "" ) then
       do radiometer = 1, database%noRadiometers
          matches(radiometer) = matches(radiometer).AND.&
            & (trim(suffix) == trim(database%radiometerInfo(radiometer)%suffix))
       end do
    end if

  end subroutine ParseRadiometerRequest

  ! -------------------------------------------  ParseBandRequest  -----

  ! This function is similar to the above one, except that more can be omitted
  ! (e.g. the spectrometer type etc.) Also the user can specify a request for
  ! upper/lower sideband.  This is reflected in the upper/lower value
  ! -1=lower, 1=upper, 0=folded.

  subroutine ParseBandRequest ( request, matches, upperLower )

    ! Dummy arguements
    character (len=*), intent(in) :: request
    logical, dimension(:) :: matches ! Result (database%noBands)
    integer, intent(out), optional :: upperLower

    ! Local variables
    character (len=len(request)) :: prefix, suffix
    character (len=1) :: thisChar, spectrometerFamilyRequested
    integer :: upperLowerRequested
    integer :: pos, bandNumberRequested, prefixLen, band
    logical :: bandStarRequested

    ! Executable code

    call SplitWords ( request, prefix, suffix, delimiter=":" )
    prefix = Capitalize(prefix)
    suffix = Capitalize(suffix)
    thisChar = prefix(1:1)
    if ( thisChar /= "B" ) call MLSMessage(MLSMSG_Error,ModuleName, &
         & "B expected in band specifier")

    ! Start with a clean slate, the user is not fussy

    upperLowerRequested = 0
    spectrometerFamilyRequested = "*"
    bandNumberRequested = 0
    bandStarRequested = .FALSE.

    prefixLen = len_trim(prefix)
    if ( prefixLen>1 ) then

       ! OK so it's not just B: something, it's more complicated

       ! Look at the second character, it's either a number or a *
       thisChar = prefix(2:2)
       if ( thisChar == "*" ) bandStarRequested = .TRUE.

       ! Now loop in from the back end to take the prefix apart
       pos = prefixLen
       ParseBandRequestParse: do
          thisChar = prefix(pos:pos)
          if ( LGE(thisChar,"0").AND.(LLE(thisChar,"9")) ) &
               & exit ParseBandRequestParse
          select case (thisChar)
          case ("*")
             if ( bandStarRequested.AND.(pos==2) ) exit ParseBandRequestParse
             ! Otherwise it's a spectrometer family request, which is already
             ! set to * so there's nothing to do.
          case ("U")
             upperLowerRequested = 1
          case ("L")
             upperLowerRequested = -1
          case default
             spectrometerFamilyRequested = thisChar
          end select
          pos = pos-1
          if (pos == 1) EXIT ParseBandRequestParse
       end do ParseBandRequestParse
       if ( pos==1 ) CALL MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Bad band request: "//TRIM(request) )

       ! Now read the band number from prefix

       if (.NOT. bandStarRequested) read (UNIT=prefix(2:pos),FMT=*) &
            & bandNumberRequested

    end if                       ! Not just B:...

    ! Now we have bandNumberRequested or 0 for don't care,
    ! UpperLower selected and spectrometerFamilyRequest
    ! Find all the bands that match our criteria.

    if ( size(matches) /= database%noBands ) call MLSMessage ( &
         & MLSMSG_Error, ModuleName, "Result wrong size" )

    matches= &
         & ( (bandNumberRequested == database%bandInfo%number) .OR. &
         &   (bandNumberRequested == 0) ) .AND. &
         & ( (spectrometerFamilyRequested == &
         &       database%bandInfo%spectrometerFamily) .OR. &
         &   (spectrometerFamilyRequested == "*") )

    ! For the suffix we have to be a little more old fasioned because we can't
    ! do a `ragged' trim

    if ( suffix /= "") then
       do band = 1, database%noBands
          matches(band) = matches(band) .AND. &
               & (TRIM(suffix) == TRIM(database%bandInfo(band)%suffix))
       end do
    end if

    if (present(upperLower)) upperLower = upperLowerRequested

  end subroutine ParseBandRequest

  ! -----------------------------------------  ParseSwitchRequest  -----

  ! This routine parses a switch request

  subroutine ParseSwitchRequest ( request, matches )

    ! Dummy arguments
    character (len=*), intent(in) :: request
    logical, dimension(:) :: matches ! (database%noSwitches)

    ! Local variables
    integer :: switch
    character (len=len(request)) :: capRequest

    ! Executable code

    capRequest=Capitalize(request)
    if (capRequest(1:1) /= "S") call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "S expected for switch spec")
    if (size(matches) /= database%noSwitches) call MLSMessage( &
         MLSMSG_Error,ModuleName,"Result wrong size")

    do switch = 1, database%noSwitches
       matches(switch) = (capRequest==database%switches(switch))
    end do
  end subroutine ParseSwitchRequest

  ! -----------------------------------  ParseSpectrometerRequest  -----

  ! This routine parses a spectrometer request. There's no need for wildcards
  ! here that make any sense.

  subroutine ParseSpectrometerRequest(request,matches)

    ! Dummy arguments
    character (len=*), intent(in) :: request
    logical, dimension(:) :: matches ! (database%noSpectrometers)

    ! Local variables
    character (len=len(request)) :: capRequest
    integer :: spectrometer

    ! Executable code

    if ( size(matches) /= database%noSpectrometers ) call MLSMessage( &
         & MLSMSG_Error, ModuleName, "Result wrong size" )

    capRequest = Capitalize(request)
    do spectrometer = 1, database%noSpectrometers
       matches(spectrometer) = (TRIM(capRequest) == &
            & TRIM(database%spectrometerInfo(spectrometer)%name))
    end do

  end subroutine ParseSpectrometerRequest

  ! ----------------------------------------  ParseChannelRequest  -----

  ! This function parses a channel request.  The form is
  ! C[<spec>+<spec>+<spec>] etc. where <spec> is n or m:n

  subroutine ParseChannelRequest ( request,firstChannel,lastChannel, &
    & channelIncluded )

    ! Dummy arguments
    character (len=*), intent(in) :: request
    integer, intent(in) :: firstChannel, lastChannel
    logical, dimension(:), pointer :: channelIncluded

    ! Local variables
    character (len=len(request)) :: field, remainder, newRemainder, &
      & firstStr, lastStr
    integer :: highestChannel, firstInRange, lastInRange
    integer :: requestLen
    character (len=1) :: keyChar ! Needed due to absoft bug

    ! Executable code

    highestChannel = firstChannel-1
    keyChar = request(1:1)
    if ( (keyChar/="C").AND.(keyChar/="c") ) call MLSMessage ( MLSMSG_Error, &
         & ModuleName, "C expected in channel specifier" )
    keyChar = request(2:2)
    if ( keyChar/="[" ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "[ expected in channel specifier" )
    requestLen = len_trim(request)
    keyChar = request(requestLen:requestLen)
    if (keyChar/="]") call MLSMessage(MLSMSG_Error,ModuleName, &
         & "] expected in channel specifier")

    remainder = Capitalize(request(3:requestLen-1))
    call allocate_test ( channelIncluded, lastChannel, "channelIncluded", &
      & ModuleName, lowBound = firstChannel, fill=.false. )

    ParseChannelParseLoop: do
       call SplitWords ( remainder, field, newRemainder, delimiter="+" )
       remainder = newRemainder
       if ( field == "" ) exit ParseChannelParseLoop

       ! See if this is a range or just a single channel
       if (index(field,":")/=0) then
          ! It's a range
          call SplitWords ( field, firstStr, lastStr, delimiter=":" )

          if ( firstStr == "" ) then
             firstInRange = firstChannel
          else
             read (UNIT=firstStr,FMT=*) firstInRange
          end if

          IF ( lastStr == "" ) then
             lastInRange = lastChannel
          else
             read (UNIT=lastStr,FMT=*) lastInRange
          endif
       else
          ! It's just a single channel
          read (UNIT=field,FMT=*) firstInRange
          lastInRange = firstInRange
       end if

       if ( firstInRange<=highestChannel ) call MLSMessage( MLSMSG_Error, &
            & ModuleName, "Channel specifier out of order" )
       channelIncluded(firstInRange:lastInRange) = .TRUE.
       highestChannel = lastInRange
    end do ParseChannelParseLoop

  end subroutine ParseChannelRequest

  ! ---------------------------------  TurnMLSChannelInfoIntoCopy  -----

  ! This routine makes sure that the channel information in a signal is a copy
  ! of the original, and doesn't merely point to it.

  subroutine TurnMLSChannelInfoIntoCopy(signals)

    ! Dummy argument
    type (MLSSignal_T), dimension(:), intent(inout) :: signals

    ! Local variables
    integer :: signal
    real(r8), dimension(:), pointer :: tempPosition, tempWidth
    logical, dimension(:), pointer :: tempIncluded

    do signal = 1, SIZE(signals)
       tempPosition => signals(signal)%channelPosition
       tempWidth => signals(signal)%channelWidth
       tempIncluded => signals(signal)%channelIncluded

       ! Don't clobber the existing one, in case it's a pointer to
       ! another signal.  Cross your fingers and hope it is; otherwise,
       ! this routine causes a memory leak.
       nullify ( signals(signal)%channelPosition, &
               & signals(signal)%channelWidth, &
               & signals(signal)%channelIncluded )

       call allocate_test ( signals(signal)%channelPosition, &
         & signals(signal)%lastChannelInBand, "channelPosition", ModuleName, &
         & lowBound = signals(signal)%firstChannelInBand )
       call allocate_test ( signals(signal)%channelWidth, &
         & signals(signal)%lastChannelInBand, "channelWidth", ModuleName, &
         & lowBound = signals(signal)%firstChannelInBand )
       call allocate_test ( signals(signal)%channelIncluded, &
         & signals(signal)%lastChannelInBand, "channelIncluded", ModuleName, &
         & lowBound = signals(signal)%firstChannelInBand )

       signals(signal)%channelPosition=tempPosition
       signals(signal)%channelWidth=tempWidth
       signals(signal)%channelIncluded=tempIncluded
    end do
  end subroutine TurnMLSChannelInfoIntoCopy

  ! --------------------------------------  ParseMLSSignalRequest  -----

  ! This subroutine takes a full request and returns an index or set of
  ! indices into the valid signals data structure, along with an indication of
  ! upper/lower sideband if relevant.  Also returned is a list of channels.

  subroutine ParseMLSSignalRequest ( request, signals, noCopy )

    ! Dummy arguments
    character (len=*), intent(in) :: request
    type (MLSSignal_T), dimension(:), pointer :: signals
    logical, optional, intent(in) :: noCopy ! Don't copy data, point to it

    ! Local paramters
    character (len=*), parameter :: OutOfOrder="Signal designation out of order"

    ! Local variables

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: upperLower
    integer :: mostAdvancedField
    character (len=len(request)) :: remainder,newRemainder,field, &
         & channelString
    logical :: radiometerMatches(database%noRadiometers), &
             & bandMatches(database%noBands), &
             & switchMatches(database%noSwitches), &
             & spectrometerMatches(database%noSpectrometers), &
             & allMatch(database%noValidSignals)
    integer :: noMatches, status, signal

    logical, dimension(:), pointer :: channelIncluded
    logical :: uniqueSpectrometerFamilyRequest
    character (len=1) :: spectrometerFamily
    logical :: useNoCopy
    character (len=1) :: keyChar ! Needed due to absoft bug

    ! Executable code

    radiometerMatches = .TRUE.
    bandMatches = .TRUE.
    switchMatches = .TRUE.
    spectrometerMatches = .TRUE.

    if ( present(noCopy) ) then
       useNoCopy = noCopy
    else
       useNoCopy = .FALSE.
    end if

    mostAdvancedField = 0
    channelString = ""

    remainder = Capitalize(request)

    ! We go through and parse the string piece by piece

    ParseMLSSignalRequestWordLoop: do
       call SplitWords ( remainder, field, newRemainder, delimiter='.' )
       if ( field == "" ) exit ParseMLSSignalRequestWordLoop
       keyChar = field(1:1)
       remainder = newRemainder
       select case ( keyChar )
       case ('R')
          IF ( mostAdvancedField /= 0 ) &
               & call MLSMessage ( MLSMSG_Error, ModuleName, &
               & OutOfOrder//TRIM(request))
          call ParseRadiometerRequest ( field, radiometerMatches )
          mostAdvancedField = 1
       case ('B')
          IF ( mostAdvancedField > 2 ) &
               & call MLSMessage ( MLSMSG_Error, ModuleName, &
               & OutOfOrder//TRIM(request))
          call ParseBandRequest ( field, bandMatches, &
               & upperLower = upperLower)
          mostAdvancedField=2
       case ('S')
          if ( mostAdvancedField > 3 ) &
               & call MLSMessage ( MLSMSG_Error,ModuleName, &
               & OutOfOrder//TRIM(request) )
          call ParseSwitchRequest ( field, switchMatches )
          mostAdvancedField = 3
       case ('C')
          channelString = field
          mostAdvancedField = 5
       case default             ! Must be spectrometer or wrong
          if ( mostAdvancedField > 4 ) &
               & call MLSMessage ( MLSMSG_Error, ModuleName, &
               & OutOfOrder//TRIM(request) )
          call ParseSpectrometerRequest ( field, spectrometerMatches )
          mostAdvancedField = 4
       end select
    end do ParseMLSSignalRequestWordLoop

    ! Now we know what they've asked for look for complete matches

    do signal = 1, database%noValidSignals

       allMatch(signal) = &
            & radiometerMatches(database%validSignals(signal)%radiometerIndex) &
            &     .AND. &
            & bandMatches(database%validSignals(signal)%bandIndex) .AND. &
            & switchMatches(database%validSignals(signal)%switchIndex) .AND. &
            & spectrometerMatches(database%validSignals(signal)%spectrometerIndex)
    end do

    noMatches = COUNT(allMatch)
    if ( noMatches==0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "No matching signal:"//TRIM(request) )

    ! Now create the result, if it's alread allocated that's an error

    if ( associated(signals) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "Signals already allocated" )
    allocate ( signals(noMatches),STAT=status )
    addr = 0
    if ( status == 0 ) then
      if ( noMatches > 0 ) addr = transfer(c_loc(signals(1)), addr)
    end if
    call test_allocate ( status, ModuleName, "signals", uBounds = noMatches, &
      & elementSize = storage_size(signals) / 8, address=addr )

    signals = pack(database%validSignals,allMatch)
    signals%upperLower = upperLower
    signals%notCopy = useNoCopy

    if ( .NOT. useNoCopy ) call TurnMLSChannelInfoIntoCopy(signals)

    ! At this point we deal with the channel specifier

    if ( trim(channelString ) /= "") then
       if ( useNoCopy ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "Cannot use no Copy with specific channels" )

       ! Also can't request specific channels when spanning more than one
       ! spectrometer type.

       if ( noMatches>1 ) then
          uniqueSpectrometerFamilyRequest = .TRUE.
          spectrometerFamily = signals(1)%spectrometerFamily
          do signal = 2, noMatches
             if ( spectrometerFamily /= signals(signal)%spectrometerFamily ) &
                  & uniqueSpectrometerFamilyRequest = .FALSE.
          end do
          if ( .NOT. uniqueSpectrometerFamilyRequest ) call MLSMessage ( &
               & MLSMSG_Error, ModuleName, &
               & "Cannot give specific channels accross multiple "// &
               & "spectrometer families" )
       end if

       call ParseChannelRequest ( channelString, &
            & signals(1)%firstchannelInBand, signals(1)%lastChannelInBand, &
            & channelIncluded )
       do signal = 1, noMatches
          signals(signal)%channelIncluded = channelIncluded
          signals(signal)%noChannelsIncluded = COUNT(channelIncluded)
       end do
       call deallocate_test ( channelIncluded, "channelIncluded", ModuleName )
    end if

  end subroutine ParseMLSSignalRequest

  ! --------------------------------------  DestroyMLSSignalsInfo  -----

  ! This routine destroys a signal or signal array as created by the
  ! ParseMLSSignalRequest routine.

  subroutine DestroyMLSSignalsInfo(signals,noError)

    ! Dummy arguments
    type (MLSSignal_T), dimension(:), pointer :: Signals
    logical, intent(in), optional :: NoError

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: s, signal, status
    logical :: useNoError

    ! Executable code

    if ( present(noError) ) then
       useNoError = noError
    else
       useNoError = .FALSE.
    endif

    if ( associated(signals) ) then
       do signal = 1, SIZE(signals)
          if ( .NOT. signals(signal)%notCopy ) then
             call deallocate_test ( signals(signal)%channelIncluded, &
               & "signals(signal)%channelIncluded", ModuleName )
             call deallocate_test ( signals(signal)%channelPosition, &
               & "signals(signal)%channelPosition", ModuleName )
             call deallocate_test ( signals(signal)%channelWidth, &
               & "signals(signal)%channelWidth", ModuleName )
          end if
       end do
       s = size(signals) * storage_size(signals) / 8
       addr = 0
       if ( s > 0 ) addr = transfer(c_loc(signals(1)), addr)
       deallocate (signals, stat=status)
       call test_deallocate ( status, ModuleName, "signals", s, address=addr )
    else
       if ( .NOT. useNoError ) call MLSMessage(MLSMSG_Error,ModuleName, &
            & "This signal not allocted" )
    end if

  end subroutine DestroyMLSSignalsInfo

  ! ----------------------------------------  UnionMLSSignalsInfo  -----

  ! This subroutine combines two sets of signal lists in a set union type of
  ! operation.  Thus duplication is avoided. Note that the result must be
  ! undefined on entry. Also note, we don't consider duplication within a or b

  subroutine UnionMLSSignalsInfo ( signalsA, signalsB, signalsUnion )

    ! Dummy arguments
    type (MLSSignal_T), dimension(:), intent(in) :: signalsA
    type (MLSSignal_T), dimension(:), intent(in) :: signalsB
    type (MLSSignal_T), dimension(:), pointer    :: signalsUnion

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: sizeA,sizeB,noSignalsInUnion
    logical :: presentInA(size(signalsB))
    integer :: status
    integer :: noRepeats

    integer :: indexA,indexB,indexUnion

    ! Executable code
    ! First work out which are the unique signals

    sizeA = SIZE(signalsA)
    sizeB = SIZE(signalsB)

    presentInA = .FALSE.

    do indexA = 1, sizeA
       presentInA = presentInA .OR. (signalsA(indexA)%signalDatabaseIndex == &
            & signalsB%signalDatabaseIndex)
    end do

    noRepeats = COUNT(presentInA)
    noSignalsInUnion = sizeA+sizeB-noRepeats

    ! Now allocate the result

    if ( associated(signalsUnion) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "Result already allocated" )
    allocate ( signalsUnion(noSignalsInUnion), stat=status )
    addr = 0
    if ( status == 0 ) then
      if ( noSignalsInUnion > 0 ) addr = transfer(c_loc(signalsUnion(1)), addr)
    end if
    call test_allocate ( status, ModuleName, "signalsUnion", &
      & uBounds = noSignalsInUnion, elementSize = storage_size(signalsUnion) / 8, &
      & address=addr )

    ! Now fill up the result

    signalsUnion(1:sizeA) = signalsA
    if ( noRepeats < sizeB ) signalsUnion(sizeA+1:noSignalsInUnion) = &
         & PACK(signalsB,(.NOT. presentInA))

    ! Copy channel info over rather than point to it

    call TurnMLSChannelInfoIntoCopy ( signalsUnion )

    ! Now we need to union any channel information. We do this using a loop,
    ! but we only need to go to sizeA

    do indexUnion = 1, sizeA
       presentInA = (signalsUnion(indexUnion)%signalDatabaseIndex == &
            & signalsB%signalDatabaseIndex)
       if ( count(presentInA) == 1 ) then
          indexB = 1
          IndexBHunt: do while ( .not. presentInA(indexB) )
             indexB = indexB + 1
          end do IndexBHunt
          signalsUnion(indexUnion)%channelIncluded = &
               & signalsUnion(indexUnion)%channelIncluded .OR. &
               & signalsB(indexB)%channelIncluded
       end if
    end do

  end subroutine UnionMLSSignalsInfo

  ! ----------------------------------  ConcatenateMLSSignalsInfo  -----

  ! This subroutine uses the one above to concatenate two lists of signals,
  ! adding signalsB into signalsA. Note that this also allows signalsA to be
  ! unallocated, in which case signalsB is just copied.

  subroutine ConcatenateMLSSignalsInfo(signalsA,signalsB)

    ! Dummy arguments
    type (MLSSignal_T), dimension(:), pointer    :: SignalsA
    type (MLSSignal_T), dimension(:), intent(in) :: SignalsB

    ! Local variables
    type (MLSSignal_T), dimension(:), pointer :: tempSignals
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: status

    ! Executable code

    if ( .NOT. associated(signalsA) ) then
       allocate ( signalsA(lbound(signalsB,1):ubound(signalsB,1)),stat=status )
       addr = 0
       if ( status == 0 ) then
         if ( size(signalsA) > 0 ) addr = transfer(c_loc(signalsA(lbound(signalsB,1))), addr)
       end if
       call test_allocate ( status, ModuleName, "signalsA", &
         & lBounds = lbound(signalsB,1), uBounds = ubound(signalsB,1), &
         & elementSize = storage_size(signalsA) / 8, address=addr )
       call TurnMLSChannelInfoIntoCopy(signalsA)
    else
       nullify ( tempSignals )
       call UnionMLSSignalsInfo ( signalsA, signalsB, tempSignals )
       call DestroyMLSSignalsInfo ( signalsA )
       signalsA => tempSignals
    end if
  end subroutine ConcatenateMLSSignalsInfo

  ! ---------------------------------  IntersectionMLSSignalsInfo  -----

  ! This subroutine is like the union one above, but this time provides the
  ! intersection.

  subroutine IntersectionMLSSignalsInfo ( signalsA, signalsB, signalsIntersection )

    ! Dummy arguments
    type (MLSSignal_T), dimension(:), pointer :: signalsA
    type (MLSSignal_T), dimension(:), pointer :: signalsB
    type (MLSSignal_T), dimension(:), pointer :: signalsIntersection

    ! Local variables
    integer :: indexA, indexB, sizeA, sizeB, status
    logical :: PresentInA(size(signalsB))
    logical, dimension(:), allocatable :: channelIncluded

    ! Executable code

    sizeA = SIZE(signalsA)
    sizeB = SIZE(signalsB)

    IF ( associated(signalsIntersection) ) call MLSMessage( &
         & MLSMSG_Error,ModuleName,"Result already allocated")

    ! Unlike the union routine this one has to be a little more do loop
    ! orientated, due to the use of channels

    do indexA = 1, sizeA
       presentInA = (signalsA(indexA)%signalDatabaseIndex == &
                   & signalsB%signalDatabaseIndex)
       if ( count(presentInA) /= 0 ) then
          indexB = 1
          IndexBHunt: do while ( .not. presentInA(indexB) )
             indexB=indexB+1
          end do IndexBHunt

          ! Now look at the channels

          allocate ( channelIncluded( &
               & signalsA(indexA)%firstChannelInBand: &
               & signalsA(indexA)%lastChannelInBand),STAT=status)
          call test_allocate ( status, ModuleName, "channelIncluded" )

          channelIncluded=signalsA(indexA)%channelIncluded .AND. &
               & signalsB(indexB)%channelIncluded

          if (count(channelIncluded) /= 0) then
             call ConcatenateMLSSignalsInfo(signalsIntersection, &
                  & signalsA(indexA:indexA))
             signalsIntersection(ubound(signalsIntersection,1))% &
                  & channelIncluded = channelIncluded
          end if
          deallocate ( channelIncluded, stat=status )
          call test_deallocate ( status, ModuleName, "channelIncluded" )
       end if
    end do

  end subroutine IntersectionMLSSignalsInfo

  ! ---------------------------------------  GetFullMLSSignalName  -----

  ! This subroutine gives the full name for an MLS signal

  subroutine GetFullMLSSignalName ( signal, fullName, rbOnly, showChannel )

    ! Dummy arguments
    type (MLSSignal_T), intent(in) :: signal
    character (len=*), intent(out) :: fullName
    logical, optional, intent(in) :: rbOnly ! Only show radiometer and band
    logical, optional, intent(in) :: showChannel ! Only show radiometer and band

    ! Local variables
    logical :: useRBOnly
    logical :: useShowChannel
    logical :: previousChanIncluded,thisChanIncluded
    integer :: channel,firstChannelInRange,rangeLen
    logical :: firstEntry
    character (len=32) :: word

    ! Executable code

    fullName=trim(signal%radiometerName)//"."//trim(signal%bandName)

    if ( present(rbOnly) ) then
       useRBOnly = rbOnly
    else
       useRBOnly = .FALSE.
    end if

    if ( present(showChannel) ) then
       useShowChannel = showChannel
    else
       useShowChannel = .FALSE.
    end if

    if ( .not. useRBOnly ) then
       fullName = trim(fullName)//"."//trim(signal%switch)//"."// &
            trim(signal%spectrometerName)
    end if

    if ( useShowChannel ) then
       fullName = trim(fullName)//".C["
       firstEntry = .TRUE.
       channel = signal%firstChannelInBand-1
       previousChanIncluded = .FALSE.
       ShowChannelLoop: do
          channel = channel+1
          if ( channel<=signal%noChannelsInBand ) then
             thisChanIncluded = signal%channelIncluded(channel)
          else
             thisChanIncluded = .FALSE.
          end if

          if ( thisChanIncluded.NEQV.previousChanIncluded ) then
             ! We have a transition of some kind
             if ( thisChanIncluded ) then ! Begining of a new range
                firstChannelInRange = channel
                write (unit=word,fmt='(i0)') firstChannelInRange
                if ( .not. firstEntry ) fullName=trim(fullName)//"+"
                firstEntry = .FALSE.
                fullName = trim(fullName)//word
             else               ! The end of a range or of array
                rangeLen = channel-firstChannelInRange
                if ( rangeLen>1 ) then
                   write (unit=word,fmt='(i0)') channel-1
                   if (rangeLen==2) then
                      fullName = trim(fullName)//"+"//word
                   else
                      fullName = trim(fullName)//":"//word
                   end if
                end if
             end if
          end if
          if ( channel==signal%lastChannelInBand+1 ) exit ShowChannelLoop
          previousChanIncluded = thisChanIncluded
       end do ShowChannelLoop
       fullName = trim(fullName)//"]"
    end if

    ! We might put channel stuff here later

  end subroutine GetFullMLSSignalName

  ! ----------------------------------------  ReadSignalsDatabase  -----

  ! This routine reads the signals database file and fills the data structure
  ! The routine doesn't do very much error checking as this file is assumed not
  ! to change very often.

  subroutine ReadSignalsDatabase ( unit )

    use Intrinsic, only: L_None
    use MLSStringLists, only: GetUniqueStrings

    ! Arguments and result

    integer, intent(in) :: unit

    ! Local parameters ------------------------------

    integer, parameter :: SDBLineLen=132
    character (len=*), parameter :: EOFMessage= &
         & "Unexpected EOF on signals database file"

    ! Some low level variables ----------------------

    character (len=SDBLineLen) :: line, first, last, rest
    integer(c_intptr_t) :: Addr         ! For tracing
    logical :: eof
    integer :: no               ! Temporary array index
    integer :: signal           ! Loop counters
    integer :: radiometer, band, switch, spectrometer, channel ! Loop counters
    integer :: firstChannel, lastChannel, noChannels
    integer :: spectrometerFamily ! Another loop counter
    integer :: wordLen          ! Length of word
    integer :: hasModifier      ! Flag
    integer :: s                ! Size in bytes of an object to deallocate
    integer :: status           ! From allocate
    integer :: index            ! General array index

    integer :: evenNo           ! For `even' channels
    real(r8) :: evenStart,evenSpacing,evenWidth ! For `even' channels

    type (MLSSignal_T), dimension(:), pointer :: tempSignal=>null()

    ! These variables are intermediate arrays to allow our database to `grow'

    type (SDBSpectrometerFamilyInfo_T) :: tempSpectrometerFamilyInfo
    character (len=SDBLineLen), dimension(:), pointer :: &
         validSignalNames

    ! These strings are a breakdown of the valid signals strings

    character (len=NameLen), dimension(:), allocatable :: validRadiometer, &
         & validBand, validSwitch, validSpectrometer
    character (len=NameLen), dimension(:), allocatable :: radiometerNames, &
         & bandNames, switchNames, spectrometerNames

    character (len=1), dimension(:), allocatable :: &
         & spectrometerFamilyChars

    integer, parameter :: MLSInstrumentNoModules=2

    character (len=3),  dimension(MLSInstrumentNoModules) :: &
       & MLSInstrumentModuleNames= (/ &
       & "GHz", &
       & "THz"/)

    ! Executable code ----------------------------------------------------
    nullify ( validSignalNames )

    ! The first section in the signals file describes the various types
    ! of spectrometer there can be.  First, we read a line and expect it to be
    ! `spectrometers'

    MLSInstrumentModuleNames(1)='GHz'
    MLSInstrumentModuleNames(2)='THz'
    !for God's sake!

    call ReadCompleteLineWithoutComments(unit,line,eof=eof)
    if (eof.OR.(Capitalize(TRIM(line))/="SPECTROMETERS")) &
         & call MLSMessage(MLSMSG_Error,ModuleName, &
         & "SPECTROMETERS expected in signals database file")

    ! Now we go through and read the family descriptions

    database%noSpectrometerFamilies=0
    SpectrometerFamilyLoop: DO
       ! Read the line
       call ReadCompleteLineWithoutComments ( unit, line, eof=eof )
       if ( eof ) call MLSMessage ( MLSMSG_Error, ModuleName, EOFMessage )
       if ( Capitalize(line)=="END" ) exit SpectrometerFamilyLoop

       ! Now parse the line

       call SplitWords ( line, first, rest, last, threeWay=.TRUE., delimiter=" " )
       tempSpectrometerFamilyInfo%name = first
       read (UNIT=rest,FMT=*) tempSpectrometerFamilyInfo%firstChannel
       read (UNIT=last,FMT=*) tempSpectrometerFamilyInfo%lastChannel
       tempSpectrometerFamilyInfo%noChannels =  &
            & tempSpectrometerFamilyInfo%lastChannel - &
            & tempSpectrometerFamilyInfo%firstChannel + 1

       database%noSpectrometerFamilies = addSpectrometerInfoToDatabase ( &
         & database%spectrometerFamilyInfo, tempSpectrometerFamilyInfo )

    end do SpectrometerFamilyLoop

    ! The next section in the file is the list of all the valid signals
    ! We'll just read these in for the moment, and parse them in the next
    ! section

    call ReadCompleteLineWithoutComments ( unit, line, eof=eof )
    if ( eof ) call MLSMessage ( MLSMSG_Error, ModuleName, EOFMessage )
    if ( Capitalize(line)/="VALID SIGNALS" ) call MLSMessage ( &
         & MLSMSG_Error, ModuleName, "Valid signals expected" )

    database%noValidSignals = 0
    ValidSignalsReadingLoop: DO
       ! Read the line
       call ReadCompleteLineWithoutComments ( unit, line, eof=eof )
       if ( eof ) call MLSMessage ( MLSMSG_Error, ModuleName, EOFMessage )
       if ( Capitalize(line)=="END" ) exit ValidSignalsReadingLoop

       database%noValidSignals = AddValidSignalNamesToDatabase ( &
         & validSignalNames, line )

    end do ValidSignalsReadingLoop

    ! The remainder of the file talks about lo frequencies etc.  We'll handle
    ! all that stuff in a minute.  First we'll take apart that information we
    ! got

    allocate ( validRadiometer(database%noValidSignals), STAT=status )
    call test_allocate ( status, ModuleName, "validRadiometer", &
      & uBounds = database%noValidSignals, &
      & elementSize = storage_size(validRadiometer) / 8 )
    allocate ( validBand(database%noValidSignals), STAT=status )
    call test_allocate ( status, ModuleName, "validBand", &
      & uBounds = database%noValidSignals, &
      & elementSize = storage_size(validBand) / 8 )
    allocate ( validSwitch(database%noValidSignals), STAT=status )
    call test_allocate ( status, ModuleName, "validSwitch", &
      & uBounds = database%noValidSignals, &
      & elementSize = storage_size(validSwitch) / 8 )
    allocate ( validSpectrometer(database%noValidSignals), STAT=status )
    call test_allocate ( status, ModuleName, "validSpectrometer", &
      & uBounds = database%noValidSignals, &
      & elementSize = storage_size(validSpectrometer) / 8 )

    allocate ( radiometerNames(database%noValidSignals), STAT=status )
    call test_allocate ( status, ModuleName, "radiometerNames", &
      & uBounds = database%noValidSignals, &
      & elementSize = storage_size(radiometerNames) / 8 )
    allocate ( bandNames(database%noValidSignals), STAT=status )
    call test_allocate ( status, ModuleName, "bandNames", &
      & uBounds = database%noValidSignals, &
      & elementSize = storage_size(bandNames) / 8 )
    allocate ( switchNames(database%noValidSignals), STAT=status )
    call test_allocate ( status, ModuleName, "switchNames", &
      & uBounds = database%noValidSignals, &
      & elementSize = storage_size(switchNames) / 8 )
    allocate ( spectrometerNames(database%noValidSignals), STAT=status )
    call test_allocate ( status, ModuleName, "spectrometerNames", &
      & uBounds = database%noValidSignals, &
      & elementSize = storage_size(spectrometerNames) / 8 )

    do signal = 1, database%noValidSignals
       ! First split into radiometer,rest,spectrometer
       call SplitWords ( validSignalNames(signal), &
            & validRadiometer(signal), rest, validSpectrometer(signal), &
            & delimiter='.',threeWay=.TRUE. )
       ! Now split rest into band,switch
       call SplitWords ( rest, validBand(signal), validSwitch(signal), &
            & delimiter='.')
    end do

    ! Now we find the unique ones of each of these and enter them in our
    ! database.

    call GetUniqueStrings ( validRadiometer, radiometerNames, &
         & database%noRadiometers )
    call GetUniqueStrings ( validBand, bandNames, database%noBands )
    call GetUniqueStrings ( validSwitch, switchNames, database%noSwitches )
    call GetUniqueStrings ( validSpectrometer, spectrometerNames, &
         & database%noSpectrometers )

    ! Now fill up our database element by element

    ! We'll go through the radiometers and fill up the radiometerInfo

    allocate ( database%radiometerInfo(database%noRadiometers), STAT=status )
    addr = 0
    if ( status == 0 ) then
      if ( database%noRadiometers > 0 ) addr = transfer(c_loc(database%radiometerInfo(1)), addr)
    end if
    call test_allocate ( status, ModuleName, "database%radiometerInfo", &
      & uBounds = database%noRadiometers, &
      & elementSize = storage_size(database%radiometerInfo) / 8, address=addr )

    do radiometer = 1, database%noRadiometers
       database%radiometerInfo(radiometer)%name = radiometerNames(radiometer)
       call SplitWords ( radiometerNames(radiometer), &
            & first, &
            & database%radiometerInfo(radiometer)%suffix, &
            & delimiter=':' )
       database%radiometerInfo(radiometer)%prefix = first

       ! We've filled the name, prefix and suffix.  Now parse the prefix into
       ! number and an optional modifer

       wordLen = LEN_TRIM(first)
       last = first(wordLen:wordLen)
       hasModifier = 0
       if ( LLT(last,"0").OR.(LGT(last,"9")) ) hasModifier=1
       read (UNIT=first(2:wordLen-hasModifier),FMT=*) &
            & database%radiometerInfo(radiometer)%number
       if ( hasModifier == 1 ) then
          database%radiometerInfo(radiometer)%modifier = last
       else
          database%radiometerInfo(radiometer)%modifier = ""
       end if
    end do

    ! Now, we similarly go through the bands and fill up that info

    allocate ( database%bandInfo(database%noBands), STAT=status )
    addr = 0
    if ( status == 0 ) then
      if ( database%noBands > 0 ) addr = transfer(c_loc(database%bandInfo(1)), addr)
    end if
    call test_allocate ( status, ModuleName, "database%bandInfo", &
      & uBounds = database%noBands, &
      & elementSize = storage_size(database%bandInfo) / 8, address=addr )

    do band = 1, database%noBands
       database%bandInfo(band)%name = bandNames(band)
       call SplitWords ( bandNames(band), &
            & first, &
            & database%bandInfo(band)%suffix, &
            & delimiter=':' )

       ! Now we parse the prefix section.  This is a B followed by a number
       ! then a spectrometer type

       wordLen = LEN_TRIM(first)
       database%bandInfo(band)%spectrometerFamily = first(wordLen:wordLen)
       read (UNIT=first(2:wordLen-1),FMT=*) &
            database%bandInfo(band)%number
    end do ! We'll sort out the spectrometer family index later

    ! Sorting out the switch settings is easy.

    allocate ( database%switches(database%noSwitches), STAT=status )
    addr = 0
    if ( status == 0 ) then
!       if ( database%noSwitches > 0 ) addr = transfer(c_loc(database%switches(1)), addr)
    end if
    call test_allocate ( status, ModuleName, "database%switches", &
      & uBounds = database%noSwitches, &
      & elementSize = storage_size(database%switches) / 8, address=addr )
    database%switches = switchNames(1:database%noSwitches)

    ! Finally we do the spectrometers

    allocate ( database%spectrometerInfo(database%noSpectrometers), &
         &  STAT=status )
    addr = 0
    if ( status == 0 ) then
      if ( database%noSpectrometers > 0 ) addr = &
        & transfer(c_loc(database%spectrometerInfo(1)), addr)
    end if
    call test_allocate ( status, ModuleName, "database%spectrometerInfo", &
      & uBounds = database%noSpectrometers, &
      & elementSize = storage_size(database%spectrometerInfo) / 8, address=addr )
    do spectrometer = 1, database%noSpectrometers
       database%spectrometerInfo(spectrometer)%name = &
            & spectrometerNames(spectrometer)
       call SplitWords ( spectrometerNames(spectrometer), &
            & database%spectrometerInfo(spectrometer)%fullFamily, &
            rest,delimiter="-" )
       read (UNIT=rest,FMT=*) &
            database%spectrometerInfo(spectrometer)%number
       database%spectrometerInfo(spectrometer)%family = &
            spectrometerNames(spectrometer)(1:1)
    end do

    ! Now as a final tidy up, we go through and fill in the family index
    ! variables in both the spectrometerInfo and bandInfo database entries.

    allocate ( spectrometerFamilyChars(database%noSpectrometerFamilies), &
         &  STAT=status )
    addr = 0
    if ( status == 0 ) then
!       if ( database%noSpectrometerFamilies > 0 ) addr = &
!         & transfer(c_loc(spectrometerFamilyChars(1)), addr)
    end if
    call test_allocate ( status, ModuleName, "spectrometerFamilyChars", &
      & uBounds = database%noSpectrometerFamilies, &
      & elementSize = storage_size(spectrometerFamilyChars) / 8, address=addr )

    do spectrometerFamily = 1, database%noSpectrometerFamilies
       spectrometerFamilyChars(spectrometerFamily) = &
            & database%spectrometerFamilyInfo(spectrometerFamily)%name(1:1)
    end do

    do spectrometer = 1, database%noSpectrometers
       database%spectrometerInfo(spectrometer)%familyIndex = &
            & LinearSearchStringArray(database%spectrometerFamilyInfo%name, &
            & database%spectrometerInfo(spectrometer)%fullFamily)
    end do

    do band = 1, database%noBands
       database%bandInfo(band)%spectrometerFamilyIndex = &
            & LinearSearchStringArray(spectrometerFamilyChars, &
            & database%bandInfo(band)%spectrometerFamily)
    end do

    ! That's basically the valid signals database filled up, no we go on and
    ! read the rest of the file which contains frequency information etc.

    ! The next section descusses the radiometer frequencies

    call ReadCompleteLineWithoutComments ( unit, line, eof=eof )
    if ( eof.OR.(Capitalize(line)/="RADIOMETER FREQUENCIES") ) &
         & call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "RADIOMETER FREQUENCIES expected in database" )

    ! Now this section is just a list of names followed by lo's

    database%radiometerInfo%lo=0.0D0

    RadiometerFrequencyLoop: DO
       call ReadCompleteLineWithoutComments ( unit, line, eof=eof )
       if ( eof ) call MLSMessage ( MLSMSG_Error, ModuleName, EOFMessage )
       if ( Capitalize(line)=="END" ) exit RadiometerFrequencyLoop
       call SplitWords ( line, first, rest, last, delimiter=' ', threeWay=.TRUE. )
       radiometer = LinearSearchStringArray ( database%radiometerInfo%name, &
            & first, caseInsensitive=.TRUE. )
       if ( radiometer == 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "No such radiometer: "//TRIM(first) )
       read (unit=rest,FMT=*) database%radiometerInfo(radiometer)%lo
       database%radiometerInfo(radiometer)%instrumentModule = &
!           & LinearSearchStringArray(MLSInstrumentModuleNames(:),last,&
            & LinearSearchStringArray(MLSInstrumentModuleNames,last,&
            & caseInsensitive=.TRUE.)

       if ( database%radiometerInfo(radiometer)%instrumentModule == &
            & L_None ) call MLSMessage ( MLSMSG_Error, &
            & ModuleName, "Unrecognised instrument module: "//last )
    end do RadiometerFrequencyLoop

    if ( MINVAL(database%radiometerInfo%lo) <= 0.0D0 ) call MLSMessage ( &
         & MLSMSG_Error, ModuleName, "Not all radiometer LOs assigned" )

    ! Now the next section is a set of similar information for the bands

    call ReadCompleteLineWithoutComments( unit, line, eof=eof )
    if ( eof.OR.(Capitalize(line)/="BAND FREQUENCIES") ) &
         & call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "BAND FREQUENCIES expected in database" )

    ! Now this section is just a list of names followed by center frequencies

    database%bandInfo%centerFreqIF = 0.0D0

    BandFrequencyLoop: do
       call ReadCompleteLineWithoutComments ( unit, line, eof=eof )
       if ( eof ) call MLSMessage ( MLSMSG_Error, ModuleName, EOFMessage )
       if (Capitalize(line)=="END" ) exit BandFrequencyLoop
       call SplitWords ( line, first, rest, delimiter=" " )
       band = LinearSearchStringArray(database%bandInfo%name, &
            & first, caseInsensitive=.TRUE.)
       if ( band == 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "No such band: "//first )
       read (unit=rest,FMT=*) database%bandInfo(band)%centerFreqIF
    end do BandFrequencyLoop

    ! Don't check here for all allocated as the WF type filters are listed
    ! seperately.

    ! Now we have the section giving the frequencies and widths of all the
    ! spectrometer channels.

    call ReadCompleteLineWithoutComments ( unit, line, eof=eof )
    if ( eof.OR.(Capitalize(TRIM(line))/="SPECTROMETER FREQUENCIES") ) &
         & call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "SPECTROMETER FREQUENCIES expected in database" )

    SpectrometerFrequencyLoop: do
       call ReadCompleteLineWithoutComments ( unit, line, eof=eof )
       if ( eof ) call MLSMessage ( MLSMSG_Error, ModuleName, EOFMessage )
       if (Capitalize(line)=="END" ) exit SpectrometerFrequencyLoop

       ! Each spectrometer family has a line describing it with a set of
       ! instructions as to how to assign the channels

       call SplitWords ( line, first, rest, delimiter=" " )
       index = LinearSearchStringArray( &
            & database%spectrometerFamilyInfo%name,first, &
            & caseInsensitive=.TRUE.)
       if ( index == 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "No such spectrometer family: "//first )

       ! Here we'll allocate the arrays for size and posiiton.  Note for the WF
       ! filters this is unneccessary, and we'll deallocate them later.
       ! But to make the code more readable I'm allocating them first.

       allocate ( database%spectrometerFamilyInfo(index)%position(&
            & database%spectrometerFamilyInfo(index)%firstChannel:&
            & database%spectrometerFamilyInfo(index)%lastChannel),&
            & STAT=status )
       addr = 0
       if ( status == 0 ) then
         if ( size(database%spectrometerFamilyInfo(index)%position) > 0 ) addr = &
           & transfer(c_loc(database%spectrometerFamilyInfo(index)%position( &
             & database%spectrometerFamilyInfo(index)%firstChannel)), addr)
       end if
       call test_allocate ( status, ModuleName, "position", &
         & lBounds = database%spectrometerFamilyInfo(index)%firstChannel, &
         & uBounds = database%spectrometerFamilyInfo(index)%lastChannel, &
         & elementSize = storage_size(database%spectrometerFamilyInfo(index)%position) / 8, &
         & address=addr )

       allocate ( database%spectrometerFamilyInfo(index)%width(&
            & database%spectrometerFamilyInfo(index)%firstChannel:&
            & database%spectrometerFamilyInfo(index)%lastChannel),&
            & STAT=status )
       addr = 0
       if ( status == 0 ) then
         if ( size(database%spectrometerFamilyInfo(index)%width) > 0 ) addr = &
           & transfer(c_loc(database%spectrometerFamilyInfo(index)%width( &
             & database%spectrometerFamilyInfo(index)%firstChannel)), addr)
       end if
       call test_allocate ( status, ModuleName, "width", &
         & lBounds = database%spectrometerFamilyInfo(index)%firstChannel, &
         & uBounds = database%spectrometerFamilyInfo(index)%lastChannel, &
         & elementSize = storage_size(database%spectrometerFamilyInfo(index)%width) / 8, &
         & address=addr )

       database%spectrometerFamilyInfo(index)%individual=.FALSE.

       select case (Capitalize(rest))
       case ('LIST') ! Positions and widths given on next lines
          read (UNIT=unit, FMT=*) &
               & database%spectrometerFamilyInfo(index)%position
          read (UNIT=unit, FMT=*) &
               & database%spectrometerFamilyInfo(index)%width

       case ('EVEN') ! No chans, start Freq, spacing and widths
          read (UNIT=unit, FMT=*) evenNo,evenStart,evenSpacing,evenWidth
          if ( evenNo /= database%spectrometerFamilyInfo(index)%noChannels ) &
               & call MLSMessage ( MLSMSG_Error, ModuleName, &
               & "Wrong number of channels for "//first )
          ! Loop and fill information
          database%spectrometerFamilyInfo(index)%width = evenWidth
          do evenNo = database%spectrometerFamilyInfo(index)%firstChannel, &
               & database%spectrometerFamilyInfo(index)%lastChannel
             database%spectrometerFamilyInfo(index)%position(evenNo) = &
                  & evenStart+evenNo*evenSpacing
          end do

       case ('INDIVIDUAL')
          database%spectrometerFamilyInfo(index)%individual=.TRUE.
          s = size(database%spectrometerFamilyInfo(index)%position) * &
            & storage_size(database%spectrometerFamilyInfo(index)%position) / 8
          addr = 0
          if ( s > 0 ) addr = transfer(c_loc( &
            & database%spectrometerFamilyInfo(index)%position(database%spectrometerFamilyInfo(index)%firstChannel)), addr)
          deallocate ( database%spectrometerFamilyInfo(index)%position, stat=status )
          call test_deallocate ( status, ModuleName, "position", s, address=addr )
          s = size(database%spectrometerFamilyInfo(index)%width) * &
            & storage_size(database%spectrometerFamilyInfo(index)%width) / 8
          addr = 0
          if ( s > 0 ) addr = transfer(c_loc( &
            & database%spectrometerFamilyInfo(index)%width( &
             & database%spectrometerFamilyInfo(index)%firstChannel)), addr)
          deallocate ( database%spectrometerFamilyInfo(index)%width, stat=status )
          call test_deallocate ( status, ModuleName, "width", s, address=addr )

       case default
          call MLSMessage ( MLSMSG_Error, ModuleName, &
               & "Unrecognized spectrometer family: "//first )
       end select
    end do SpectrometerFrequencyLoop

    ! The final section of the file is the section dealing with specific
    ! channels. Here, having built up most of the database, we can find the
    ! relevant data using our standard parse routines.

    ! Now we have read all the gory information from our database
    ! we need to fill up our valid signals database

    allocate ( database%validSignals(database%noValidSignals),STAT=status )
    addr = 0
    if ( status == 0 ) then
      if ( database%noValidSignals > 0 ) addr = &
        & transfer(c_loc(database%validSignals(1)), addr)
    end if
    call test_allocate ( status, ModuleName, "database%validSignals", &
      & uBounds = database%noValidSignals, &
      & elementSize = storage_size(database%validSignals) / 8, address=addr )

    do signal=1,database%noValidSignals

       database%validSignals(signal)%signalDatabaseIndex=signal

       ! Deal with radiometer info first.

       radiometer = LinearSearchStringArray( &
            & database%radiometerInfo%name,validRadiometer(signal))
       database%validSignals(signal)%radiometerIndex = radiometer
       database%validSignals(signal)%instrumentModule = &
            & database%radiometerInfo(radiometer)%instrumentModule
       database%validSignals(signal)%radiometerName = &
            & database%radiometerInfo(radiometer)%name
       database%validSignals(signal)%radiometerPrefix = &
            & database%radiometerInfo(radiometer)%prefix
       database%validSignals(signal)%radiometerSuffix = &
            & database%radiometerInfo(radiometer)%suffix
       database%validSignals(signal)%radiometerNumber = &
            & database%radiometerInfo(radiometer)%number
       database%validSignals(signal)%radiometerModifier = &
            & database%radiometerInfo(radiometer)%modifier
       database%validSignals(signal)%lo = &
            & database%radiometerInfo(radiometer)%lo

       ! Now the band info

       band = LinearSearchStringArray(database%bandInfo%name,validBand(signal))
       database%validSignals(signal)%bandIndex = band
       database%validSignals(signal)%bandName = database%bandInfo(band)%name
       database%validSignals(signal)%bandSuffix = database%bandInfo(band)%suffix
       database%validSignals(signal)%spectrometerFamily = &
            & database%bandInfo(band)%spectrometerFamily
       database%validSignals(signal)%bandCenterFreqIF = &
            & database%bandInfo(band)%centerFreqIF

       ! Now the switch info

       switch = LinearSearchStringArray(database%switches,validSwitch(signal))
       database%validSignals(signal)%switchIndex = switch
       database%validSignals(signal)%switch = database%switches(switch)

       ! Now the spectrometer info

       spectrometer = LinearSearchStringArray(database%spectrometerInfo%name, &
            & validSpectrometer(signal))
       spectrometerFamily = database%spectrometerInfo(spectrometer)%familyIndex

       database%validSignals(signal)%spectrometerIndex = spectrometer
       database%validSignals(signal)%spectrometerFamilyIndex = spectrometerFamily

       database%validSignals(signal)%spectrometerName = &
            & database%spectrometerInfo(spectrometer)%name
       database%validSignals(signal)%fullSpectrometerFamily = &
            & database%spectrometerInfo(spectrometer)%fullFamily
       database%validSignals(signal)%spectrometerNumber = &
            & database%spectrometerInfo(spectrometer)%number

       ! Now the information on the channels

       firstChannel = database% &
            & spectrometerFamilyInfo(spectrometerFamily)%firstChannel
       lastChannel = database% &
            & spectrometerFamilyInfo(spectrometerFamily)%lastChannel
       noChannels = database% &
            & spectrometerFamilyInfo(spectrometerFamily)%noChannels

       ! Now we compute the channel positions in IF space

       database%validSignals(signal)%firstChannelInBand = firstChannel
       database%validSignals(signal)%lastChannelInBand = lastChannel
       database%validSignals(signal)%noChannelsInBand = noChannels

       allocate ( database%validSignals(signal)%channelPosition &
            & (firstChannel:lastChannel),STAT=status )
       addr = 0
       if ( status == 0 ) then
         if ( size(database%validSignals(signal)%channelPosition) > 0 ) addr = &
           & transfer(c_loc(database%validSignals(signal)%channelPosition(firstChannel)), addr)
       end if
       call test_allocate ( status, ModuleName, "channelPosition", &
         & lBounds = firstChannel, uBounds = lastChannel, &
         & elementSize = storage_size(database%validSignals(signal)%channelPosition) / 8, &
         & address=addr )

       allocate ( database%validSignals(signal)%channelWidth &
            & (firstChannel:lastChannel),STAT=status )
       addr = 0
       if ( status == 0 ) then
         if ( size(database%validSignals(signal)%channelWidth) > 0 ) addr = &
           & transfer(c_loc(database%validSignals(signal)%channelWidth(firstChannel)), addr)
       end if
       call test_allocate ( status, ModuleName, "channelWidth", &
         & lBounds = firstChannel, uBounds = lastChannel, &
         & elementSize = storage_size(database%validSignals(signal)%channelWidth) / 8, &
         & address=addr )

       allocate ( database%validSignals(signal)%channelIncluded &
            & (firstChannel:lastChannel),STAT=status )
       addr = 0
       if ( status == 0 ) then
         if ( size(database%validSignals(signal)%channelIncluded) > 0 ) addr = &
           & transfer(c_loc(database%validSignals(signal)%channelIncluded(firstChannel)), addr)
       end if
       call test_allocate ( status, ModuleName, "channelIncluded", &
         & lBounds = firstChannel, uBounds = lastChannel, &
         & elementSize = storage_size(database%validSignals(signal)%channelIncluded) / 8, &
         & address=addr )

       database%validSignals(signal)%channelIncluded = .TRUE.

       do channel = firstChannel, lastChannel
          ! We'll have to defer dealing with the individual channels
          if ( .NOT. database%spectrometerFamilyInfo(spectrometerFamily)% &
               & individual ) then
             database%validSignals(signal)%channelPosition(channel) = &
                  & database%validSignals(signal)%bandCenterFreqIF + &
                  & database%spectrometerFamilyInfo(spectrometerFamily)% &
                  & position(channel)
             database%validSignals(signal)%channelWidth(channel) = &
                  & database%spectrometerFamilyInfo(spectrometerFamily)% &
                  & width(channel)
          end if
       end do                   ! Channel loop
    end do                      ! Signal loop

    ! The next section of the file deals with the WF type channels.  This is
    ! actually very easy.  We simply read the spec and parse it as we've done
    ! others and store it in our database, making use of the
    ! ParseMLSSignalRequest routine which now works because we've fill up the
    ! daabase.

    call ReadCompleteLineWithoutComments ( unit, line, eof=eof )
    if ( eof.OR.(Capitalize(line)/="SPECIFIC CHANNELS") ) &
         & call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "SPECIFIC CHANNELS expected in database" )

    SpecificChannelLoop: do
       call ReadCompleteLineWithoutComments ( unit, line, eof=eof )
       if ( eof ) call MLSMessage ( MLSMSG_Error, ModuleName, EOFMessage )

       ! This line will either be END or a spec followed by the word 'List'
       line = Capitalize(line)
       if ( line=="END" ) exit SpecificChannelLoop

       call SplitWords ( line, first, last, delimiter=" " )
       if ( last/="LIST") call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "List expected for specific channel information" )

       ! This will be followed by a list of positions and widths, we'll read
       ! this directory into our database by using the noCopy option.

       call ParseMLSSignalRequest ( first, tempSignal, noCopy=.TRUE. )
       if ( SIZE(tempSignal)/=1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Ambiguous spec.:"//first )

       read (UNIT=unit,FMT=*) tempSignal(1)%channelPosition, &
            & tempSignal(1)%channelWidth

       ! Now destroy this temporary information

       call DestroyMLSSignalsInfo ( tempSignal )

    end do SpecificChannelLoop

    ! Now we tidy up our arrays and exit

    s = size(validSignalNames) * storage_size(validSignalNames) / 8
    addr = 0
!     if ( s > 0 ) addr = transfer(c_loc(validSignalNames(1)), addr)
    deallocate( validSignalNames, stat=status )
    call test_deallocate ( status, ModuleName, "validSignalNames", s, address=addr )
    s = size(radiometerNames) * storage_size(radiometerNames) / 8
    deallocate ( radiometerNames, stat=status )
    call test_deallocate ( status, ModuleName, "radiometerNames", s )
    s = size(bandNames) * storage_size(bandNames) / 8
    deallocate ( bandNames, stat=status )
    call test_deallocate ( status, ModuleName, "bandNames", s )
    s = size(switchNames) * storage_size(switchNames) / 8
    deallocate ( switchNames, stat=status )
    call test_deallocate ( status, ModuleName, "switchNames", s )
    s = size(spectrometerNames) * storage_size(switchNames) / 8
    deallocate ( spectrometerNames, stat=status )
    call test_deallocate ( status, ModuleName, "spectrometerNames", s )
    s = size(spectrometerFamilyChars) * storage_size(spectrometerFamilyChars) / 8
    deallocate ( spectrometerFamilyChars, stat=status )
    call test_deallocate ( status, ModuleName, "spectrometerFamilyChars", s )

  contains

    integer function AddSpectrometerInfoToDatabase ( Database, Item )
      use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

      ! Dummy arguments
      type (SDBSpectrometerFamilyInfo_T), dimension(:), pointer :: DATABASE
      type (SDBSpectrometerFamilyInfo_T), intent(in) :: ITEM

      ! Local variables
      type (SDBSpectrometerFamilyInfo_T), dimension(:), pointer :: tempDatabase

      include "addItemToDatabase.f9h" 

      addSpectrometerInfoToDatabase = newSize

    end function AddSpectrometerInfoToDatabase

    integer function AddValidSignalNamesToDatabase ( Database, Item )
      use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

      ! Dummy arguments
      character (len=SDBLineLen), dimension(:), pointer :: DATABASE
      character (len=SDBLineLen), intent(in) :: ITEM

      ! Local variables
      character (len=SDBLineLen), dimension(:), pointer :: tempDatabase

      include "addItemToDatabase.f9h" 

      addValidSignalNamesToDatabase = newSize

    end function AddValidSignalNamesToDatabase

  end subroutine ReadSignalsDatabase

  ! -------------------------------------  DestroySignalsDatabase  -----

  ! This routine deallocates all the information dealing with the signals
  ! database

  subroutine DestroySignalsDatabase

   ! Local variable

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: i, s, status

    ! Executable code

    database%noRadiometers = 0
    database%noBands = 0
    database%noSwitches = 0
    database%noSpectrometers = 0
    database%noSpectrometerFamilies = 0
    database%noValidSignals = 0

    s = size(database%radiometerInfo) * storage_size(database%radiometerInfo) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(database%radiometerInfo(1)), addr)
    deallocate ( database%radiometerInfo, stat=status )
    call test_deallocate ( status, ModuleName, "database%radiometerInfo", s, &
      & address=addr )
    s = size(database%bandInfo) * storage_size(database%bandInfo) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(database%bandInfo(1)), addr)
    deallocate ( database%bandInfo, stat=status )
    call test_deallocate ( status, ModuleName, "database%bandInfo", s, &
      & address=addr )
    s = size(database%switches) * storage_size(database%switches) / 8
!     if ( s > 0 ) addr = transfer(c_loc(database%switches(1)), addr)
    deallocate ( database%switches, stat=status )
    call test_deallocate ( status, ModuleName, "database%switches", s, &
      & address=addr )
    s = size(database%spectrometerInfo) * storage_size(database%spectrometerInfo) / 8
    if ( s > 0 ) addr = transfer(c_loc(database%spectrometerInfo(1)), addr)
    deallocate ( database%spectrometerInfo, stat=status )
    call test_deallocate ( status, ModuleName, "database%spectrometerInfo", s, &
      & address=addr )

    do i = 1, database%noSpectrometerFamilies
       s = size(database%spectrometerFamilyInfo(i)%position) * &
         & storage_size(database%spectrometerFamilyInfo(i)%position) / 8
       addr = 0
       if ( s > 0 ) addr = transfer(c_loc( &
         & database%spectrometerFamilyInfo(i)%position(database%spectrometerFamilyInfo(i)%firstChannel)), addr)
       deallocate ( database%spectrometerFamilyInfo(i)%position, stat=status )
       call test_deallocate ( status, ModuleName, &
         & "database%spectrometerFamilyInfo(i)%position", s, address=addr )
       s = size(database%spectrometerFamilyInfo(i)%width) * &
         & storage_size(database%spectrometerFamilyInfo(i)%width) / 8
       addr = 0
       if ( s > 0 ) addr = transfer(c_loc( &
         & database%spectrometerFamilyInfo(i)%width(database%spectrometerFamilyInfo(i)%firstChannel)), addr)
       deallocate ( database%spectrometerFamilyInfo(i)%width, stat=status )
       call test_deallocate ( status, ModuleName, &
         & "database%spectrometerFamilyInfo(i)%width", s, address=addr )
    end do
    s = size(database%spectrometerFamilyInfo) * &
      & storage_size(database%spectrometerFamilyInfo) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc( &
      & database%spectrometerFamilyInfo(i)%width(database%spectrometerFamilyInfo(i)%firstChannel)), addr)
    deallocate ( database%spectrometerFamilyInfo, stat=status)
    call test_deallocate ( status, ModuleName, "database%spectrometerFamilyInfo", s, &
      & address=addr )

    do i = 1, database%noValidSignals
       s = size(database%validSignals(i)%channelIncluded) * &
         & storage_size(database%validSignals(i)%channelIncluded) / 8
       addr = 0
       if ( s > 0 ) addr = transfer(c_loc( &
         & database%validSignals(i)%channelIncluded(lbound(database%validSignals(i)%channelIncluded,1))),addr)
       deallocate ( database%validSignals(i)%channelIncluded, stat=status )
       call test_deallocate ( status, ModuleName, &
         & "database%validSignals(i)%channelIncluded", s, address=addr )
       s = size(database%validSignals(i)%channelPosition) * &
         & storage_size(database%validSignals(i)%channelPosition) / 8
       addr = 0
       if ( s > 0 ) addr = transfer(c_loc( &
         & database%validSignals(i)%channelPosition(lbound(database%validSignals(i)%channelPosition,1))),addr)
       deallocate ( database%validSignals(i)%channelPosition, stat=status )
       call test_deallocate ( status, ModuleName, &
         & "database%validSignals(i)%channelPosition", s, address=addr )
       s = size(database%validSignals(i)%channelWidth) * &
         & storage_size(database%validSignals(i)%channelWidth) / 8
       addr = 0
       if ( s > 0 ) addr = transfer(c_loc( &
         & database%validSignals(i)%channelWidth(lbound(database%validSignals(i)%channelWidth,1))),addr)
       deallocate ( database%validSignals(i)%channelWidth, stat=status )
       call test_deallocate ( status, ModuleName, &
         & "database%validSignals(i)%channelWidth", s, address=addr )
    end do
    s = size(database%validSignals) * &
      & storage_size(database%validSignals) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(database%validSignals(1)), addr)
    deallocate ( database%validSignals, stat=status )
    call test_deallocate ( status, ModuleName, "database%validSignals", s, address=addr )

  end subroutine DestroySignalsDatabase

  ! --------------------------------------  GetMLSRadiometerNames  -----

  ! This routine returns an array of the MLS radiometer names from the database

  subroutine GetMLSRadiometerNames ( names )

    ! Dummy arguments
    character (len=NameLen), dimension(:), pointer :: names

    ! Local variables
    integer :: status

    ! Executable code

    call allocate_test ( names, database%noRadiometers, "MLSRadiometerNames", &
      & moduleName )

    names = database%radiometerInfo%name
  end subroutine GetMLSRadiometerNames

  ! --------------------------------------------  GetMLSBandNames  -----

  ! This routine returns an array of names of bands in MLS from the database.

  subroutine GetMLSBandNames ( names )

    ! Dummy arguments
    character (len=NameLen), dimension(:), pointer :: names

    ! Local variables
    integer :: status

    ! Executable code

    call allocate_test ( names, database%noBands, "MLSBandNames", moduleName )

    names = database%bandInfo%name
  end subroutine GetMLSBandNames

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

end module MLSSignalNomenclature
!=============================================================================

!
! $Log$
! Revision 2.16  2017/01/04 17:35:32  pwagner
! Pre-nullify pointer components of MLSSignalsDatabase_T to avoid crashes during addItem ..
!
! Revision 2.15  2016/03/08 21:08:18  pwagner
! Fixed two errors causing some crashes
!
! Revision 2.14  2015/03/28 01:18:07  vsnyder
! Some spiffing.  Some reorganization.
! Added stuff to trace allocate/deallocate addresses -- some commented out
! because NAG build 1017 doesn't yet allow arrays as arguments to C_LOC.
!
! Revision 2.13  2014/09/30 16:31:02  pwagner
! Uses Allocate_test, etc. from Allocate_Deallocate, back to module-wise scope
!
! Revision 2.12  2014/09/29 22:50:26  vsnyder
! Add ONLY clause to USEs, move some down to subprogram scope
!
! Revision 2.11  2014/09/05 21:58:20  pwagner
! Remove wrong call to allocate_test with signals array
!
! Revision 2.10  2014/09/05 00:09:45  vsnyder
! More complete and accurate allocate/deallocate size tracking.
! Convert some local pointer temps to automatic.  Some cannonball polishing.
!
! Revision 2.9  2013/06/13 21:04:58  vsnyder
! Don't look at noError without checking whether it's present
!
! Revision 2.8  2013/06/12 02:12:07  vsnyder
! Cruft removal
!
! Revision 2.7  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.6  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.5  2004/08/04 23:19:01  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.4  2002/10/08 00:09:12  pwagner
! Added idents to survive zealous Lahey optimizer
!
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
