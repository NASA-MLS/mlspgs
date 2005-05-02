! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Parse_Signal_M

  implicit NONE
  private
! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

! Expand_Signal_List      Convert subtree to fully described signals
! Get_Individual_Signals  Expand signal array to fully-qualified signals
! Parse_Signal            Return database indices matching signal string
! === (end of toc) ===

! === (start of api) ===
! Expand_Signal_List ( int node, signal_t* theseSignals[], log* error )
! Get_Individual_Signals ( char* outSignals[], int inSignals_ind[] char* inSignals_txt[] )
! Parse_Signal  ( char* Signal_String, int* Signal_Indices(:),
!    [int Tree_Index], [int Sideband],[log* Channels(:)], [int OnlyCountEm] )
! === (end of api) ===

  public :: Parse_Signal, Expand_Signal_List, Get_Individual_Signals

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains

  ! -----------------------------------------  Expand_signal_list  -----
  subroutine Expand_signal_list ( node, theseSignals, error )
    ! This uses parse_signal below to convert an array of strings
    ! at tree node 'node' into a fully described array of signals
    ! populated with channel flags and everything

    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use MLSSignals_m, only: SIGNALS, SIGNAL_T, DESTROYSIGNALDATABASE
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_ERROR
    use Tree, only: SUBTREE, SUB_ROSA, NSONS
    use String_table, only: GET_STRING

    integer, intent(in) :: NODE         ! Tree node
    type ( Signal_T ), pointer, dimension(:) :: THESESIGNALS ! Result
    logical, intent(out) :: ERROR

    ! Local variables
    integer :: J                        ! Loop counter
    integer :: STATUS                   ! Flag from allocate etc.
    character (len=132) :: SIGNALSTRING ! One signal name
    integer, dimension(:), pointer :: SIGNALINDS ! Array of signal indicies
    logical, dimension(:), pointer :: CHANNELS ! Which channels in signal
    integer :: SIDEBAND                 ! One sideband
    integer :: Son                      ! of Node
    integer :: WANTED                   ! Which of the possible matches do we want

    ! Executable code

    ! Make sure the result is deallocated, do some setup
    call DestroySignalDatabase ( theseSignals, justChannels=.true. )
    allocate ( theseSignals (nsons(node)-1), stat = status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, MLSMSG_Allocate // 'theseSignals' )
    nullify ( channels, signalInds )
    error = .false.

    ! Loop over the strings
    do j = 1, nsons(node)-1
      son = subtree(j+1,node)
      call get_string ( sub_rosa(son), signalString, strip=.true.)
      call parse_Signal ( signalString, signalInds, &
        & tree_index=son, sideband=sideband, channels=channels )
      if ( .not. associated(signalInds) ) then ! A parse error occurred
        error = .true.
        ! Could tidy up, but expect to crash anyway
        return
      end if
      ! parse_signal may have returned several signals (e.g. different
      ! switch settings). Later on choose the `right' one from the match
      ! For the moment choose the first !????
      wanted = 1
      theseSignals(j) = signals(signalInds(wanted))
      theseSignals(j)%sideband = sideband
      ! Don't hose channels in database, though shouldn't be an issue
      nullify ( theseSignals(j)%channels ) 
      
      call allocate_Test ( theseSignals(j)%channels, &
        & ubound(theseSignals(j)%frequencies,1), 'signals%channels', &
        & ModuleName, lowBound=lbound(theseSignals(j)%frequencies,1) )
      if ( associated(channels) ) then
        ! The reason that this is so messy is that channels has l/u bounds
        ! of the highest and lowest channels asked for.
        theseSignals(j)%channels(:lbound(channels,1)-1) = .false.
        theseSignals(j)%channels(lbound(channels,1):ubound(channels,1)) = &
          channels
        theseSignals(j)%channels(ubound(channels,1)+1:) = .false.
      else
        theseSignals(j)%channels = .true.
      end if
      call deallocate_test ( channels, 'channels', ModuleName )
      call deallocate_test ( signalInds, 'signalInds', ModuleName )
    end do                          ! End loop over listed signals

    ! Executable code
  end subroutine Expand_signal_list

  ! -------------------------------------  Get_Individual_Signals  -----
  subroutine Get_Individual_Signals ( OutSignals, InSignals_Ind, InSignals_Txt )
  ! Given InSignals_Ind or InSignals_Txt, create OutSignals where each of
  ! OutSignals is a signal denoted in InSignals_*, but describes exactly
  ! one radiometer, band, switch, spectrometer, sideband and channel. 
  ! OutSignals is allocated using Allocate_Test, so don't send in an
  ! undefined pointer! There will be no duplicates in OutSignals. 
  ! Zero-size InSignals_* works.

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_m, only: GetSignalName, MaxSigLen, Signals
    use Sort_m, only: Sort
    use String_Table, only: Get_String

    character(len=*), pointer :: OutSignals(:)
    integer, intent(in), optional :: InSignals_Ind(:)
    character(len=*), intent(in), optional :: InSignals_Txt(:)

    logical, pointer :: Channels(:)
    integer :: I, J, K, L, N
    character(len=len(outSignals)), pointer :: MySignals(:)
    character(len=maxSigLen) :: OneSignal
    integer :: Sideband
    integer, pointer :: SignalIndices(:)

    nullify ( channels, mySignals, signalIndices )
    ! Determine the size of mySignals
    l = 0
    if ( present(inSignals_ind) ) n = size(inSignals_ind)
    if ( present(inSignals_txt) ) n = size(inSignals_txt)
    do i = 1, n
      if ( present(inSignals_ind) ) &
        & call get_string ( inSignals_ind(i), oneSignal, strip=.true. )
      if ( present(inSignals_txt) ) &
        & oneSignal = inSignals_txt(i)
      call parse_signal ( oneSignal, signalIndices, sideband=sideband, &
        & channels=channels )
      if ( .not. associated(signalIndices) ) call MLSMessage ( MLSMSG_Error, &
        & moduleName, 'Unable to parse signal ' // trim(oneSignal) )
      if ( .not. associated(channels) ) then
        ! All channels are desired.  Assume channels available in all
        ! selected signals are the same.
        k = size(signals(signalIndices(1))%frequencies)
      else
        k = count(channels)
      end if
      l = l + k
      do j = 1, size(signalIndices)
        l = l + k
        if ( sideband == 0 .and. signals(signalIndices(j))%singleSideband == 0 ) &
          & l = l + k ! pick up both sidebands
      end do ! j = 1, size(signalIndices)
    end do ! i = 1, n
    ! Fill mySignals
    call allocate_test ( mySignals, l, 'mySignals', moduleName )
    if ( l /= 0 ) then
      l = 0
      do i = 1, n
        if ( present(inSignals_ind) ) &
          & call get_string ( inSignals_ind(i), oneSignal, strip=.true. )
        if ( present(inSignals_txt) ) &
          & oneSignal = inSignals_txt(i)
        call parse_signal ( oneSignal, signalIndices, sideband=sideband, &
          & channels=channels )
        if ( .not. associated(channels) ) then
          ! All channels are desired.  Assume channels available in all
          ! selected signals are the same.
          call allocate_test ( channels, &
            & ubound(signals(signalIndices(1))%frequencies,1), 'myChannels', &
            & moduleName, &
            & lowBound=lbound(signals(signalIndices(1))%frequencies,1) )
          channels = .true.
        end if
        do j = 1, size(signalIndices)
          do k = lbound(channels,1), ubound(channels,1)
            if ( channels(k) ) then
              l = l + 1
              if ( sideband == 0 .and. signals(signalIndices(j))%singleSideband == 0 ) then
                call getSignalName ( signalIndices(j), mySignals(l), &
                  & channel=k, sideband=-1 )
                l = l + 1
                call getSignalName ( signalIndices(j), mySignals(l), &
                  & channel=k, sideband=+1 )
              else
                call getSignalName ( signalIndices(j), mySignals(l), &
                  & channel=k, sideband=sideband )
              end if
            end if
          end do ! k = lbound(channels), ubound(channels)
        end do ! j = 1, size(signalIndices)
      end do ! i = 1, n
      ! Delete duplicate signals
      call sort ( mySignals, 1, l )
      k = 1
      do i = 2, l
        if ( mySignals(i) /= mySignals(i-1) ) k = k + 1
      end do
    else
      k = 0
    end if
    ! Copy mySignals to outSignals, deleting duplicates
    call allocate_test ( outSignals, k, 'Outsignals', moduleName )
    if ( l /= 0 ) then
      outSignals(1) = mySignals(1)
      k = 1
      do i = 2, l
        if ( mySignals(i) /= mySignals(i-1) ) then
          k = k + 1
          outSignals(k) = mySignals(i)
        end if
      end do
    end if
    ! Clean up
    call deallocate_test ( channels, 'Channels', moduleName )
    call deallocate_test ( mySignals, 'MySignals', moduleName )
    call deallocate_test ( signalIndices, 'SignalIndices', moduleName )
  end subroutine Get_Individual_Signals

  ! -----------------------------------------------  Parse_signal  -----

  subroutine Parse_Signal ( Signal_String, Signal_Indices, &
    & Tree_Index, Sideband, Channels, OnlyCountEm )

  ! Parse a signal string.  Return the indices in the signal database of
  ! signals that match the signal string.

  ! A signal string is of the form:
  ! [Radiometer[:suffix].] [Band[:suffix].] [S switch-number.] [Spectrometer.]
  ! [C channel_number[:|+channel_number]*

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Declaration_Table, only: Decls, Get_Decl, Label
    use Init_MLSSignals_m, only: S_Band, S_Radiometer, S_SpectrometerType
    use Intrinsic, only: Spec_Indices
    use Lexer_Core, only: Print_Source, Token
    use Lexer_m, only: Lex_Signal
    use MLSSignals_m, only: Bands, Radiometers, Signals
    use MoreTree, only: Get_Spec_ID
    use Output_m, only: Output
    use String_table, only: Display_string, Get_String
    use Symbol_Table, only: Dump_Symbol_Class, Enter_Terminal
    use Symbol_Types, only: T_Colon, T_Dot, T_End_of_input, T_Identifier, &
      & T_Minus, T_Plus
    use Tree, only: Decoration, Source_Ref

    character(len=*), intent(in) :: Signal_String ! Input
    integer, pointer :: Signal_Indices(:)         ! Indices in the signals
    ! database.  Deallocated here using Allocate_Test, so don't start with an
    ! undefined pointer!  Upon return, if Signal_Indices is not associated
    ! (and OnlyCountEm is not present), an error occurred.
    integer, intent(in), optional :: Tree_Index   ! To get line and column
    !                                               numbers for messages
    integer, intent(out), optional :: Sideband    ! Zero if no band is
    ! specified, or if a band is specified but it has none of [UuLl] in it.
    ! +1 if a band is specified and it has [Uu] in it.  -1 if a band is
    ! specified and it has [Ll] in it.
    logical, pointer, dimension(:), optional :: Channels    ! Output: Channel
    ! numbers mentioned after the .C part of the signal specification. They
    ! are allocated here using Allocate_Test, so don't start with an
    ! undefined pointer!
    integer, intent(out), optional :: OnlyCountEm ! If OnlyCountEm is present,
    ! Signal_Indices is NOT touched; instead, OnlyCountEm returns the number
    ! of matching signals.

    integer :: Band_i         ! Database index
    type(decls) :: Decl       ! Declaration of an identifier
    integer :: Error          ! >0 indicates an error
    integer :: I              ! Temporary, loop inductor, subscript
    integer, dimension(len(signal_string)/2,2) :: MyChannelNumbers
    logical, pointer, dimension(:) :: MyChannels
    integer :: MySideBand     ! Index of [UuLl] in band, then 0 or +/- 1.
    integer :: Next           ! One of Saw... parameters below, to indicate
    !                           the category of the next identifier
    integer :: NumChannelNums ! How many elements of MyChannelNumbers are used
    integer :: Radiometer_i   ! Database index
    integer :: SawWhat        ! What has been seen so far.  See Saw... below
    integer :: Spectrometer_i ! Database index
    integer :: Spectrometer_n ! Number after "-", in e.g. FB25-3
    integer :: Status         ! From an internal READ
    integer :: Switch         ! Switch number, from the token text
    character(len=32) :: TestText   ! From a "%suffix" field via sub_rosa
    type(token) :: The_Token  ! that results from lexing the signal string
    character(len=32) :: TokenText  ! Duh
    integer :: Where     ! The current analysis point in the signal string

    ! States.  The pieces of the radiometer specification must be presented
    ! in order, but pieces can be missing.
    integer, parameter :: SawNothing = 0
    integer, parameter :: SawRadiometer = sawNothing + 1
    integer, parameter :: SawBand = sawRadiometer + 1
    integer, parameter :: SawSwitch = sawBand + 1
    integer, parameter :: SawSpectrometerType = sawSwitch + 1
    integer, parameter :: SawChannel = sawSpectrometerType + 1
    integer, parameter :: SawChannelColon = sawChannel + 1
    integer, parameter :: SawChannelPlus = SawChannelColon + 1

    ! Parameters for error message codes
    integer, parameter :: ChannelsDecrease = 1       ! Channel numbers not in order
    integer, parameter :: Invalid = ChannelsDecrease + 1 ! Invalid combination
    integer, parameter :: Malformed = invalid + 1    ! Expected Letter Digits+
    integer, parameter :: OutOfOrder = malformed + 1 ! Identifiers not in order
    integer, parameter :: NoSuffix = outOfOrder + 1  ! Suffix not allowed
    integer, parameter :: NotEnough = noSuffix + 1   ! Not enough info
    integer, parameter :: NotLabel = notEnough + 1   ! Not the label of ...
    integer, parameter :: NotNumber = notLabel + 1   ! Expected a number
    integer, parameter :: Unexpected = notNumber + 1 ! Unexpected token
    integer, parameter :: WrongSuffix = unexpected + 1

    ! Initialize
    band_i = -1
    error = 0
    nullify ( myChannels )
    mySideband = 0
    next = sawNothing
    numChannelNums = 0
    radiometer_i = -1
    sawWhat = sawNothing
    call deallocate_test ( signal_indices, moduleName, "Signal_Indices" )
    spectrometer_i = -1
    spectrometer_n = -1
    switch = -1
    where = 0

    ! Get the tokens; verify they are the correct kinds of thing-o's.
o:  do
      call lex_signal ( signal_string, where, the_token )
      if ( the_token%class == t_end_of_input ) &
    exit o
      if ( the_token%class /= t_identifier ) &
        & call announce_error ( unexpected, signal_string, tree_index, &
          & expected = (/ t_identifier /) )
      call get_string ( the_token%string_index, tokenText, strip=.true. )
      if ( verify(tokenText(1:1), 'CcSs') == 0 ) then
        read ( tokenText(2:), *, iostat=status ) i
        if ( status /= 0 ) &
          & call announce_error ( malformed, signal_string, tree_index )
        if ( verify(tokenText(1:1), 'Ss') == 0 ) then
          next = sawSwitch
          switch = i
        else if ( verify(tokenText(1:1), 'Cc') == 0 ) then
          next = sawChannel
          numChannelNums = 1
          MyChannelNumbers(1,:) = i
          ! Process [":"<number>](("+"<number>)+[":"<number>])*
          do
            call lex_signal ( signal_string, where, the_token )
            if ( the_token%class == t_end_of_input ) &
    exit o
            if ( the_token%class == t_colon ) then
              next = sawChannelColon
            else if ( the_token%class == t_plus ) then
              next = sawChannelPlus
            else
              call announce_error ( unexpected, signal_string, tree_index, &
                & expected = (/ t_colon, t_plus /) )
            end if
            if ( next <= sawWhat ) &
              & call announce_error ( outOfOrder, signal_string, tree_index )
            call lex_signal ( signal_string, where, the_token )
            if ( the_token%class /= t_identifier ) then
              call announce_error ( unexpected, signal_string, tree_index, &
                & expected = (/ t_identifier /) )
          exit
            end if
            call get_string ( the_token%string_index, tokenText, strip=.true. )
            read ( tokenText, *, iostat=status ) i
            if ( status /= 0 ) &
              & call announce_error ( malformed, signal_string, tree_index )
            if ( next == sawChannelColon ) then
              if ( i < myChannelNumbers(numChannelNums,1) ) &
                & call announce_error ( channelsDecrease, signal_string, &
                  & tree_index )
              myChannelNumbers(numChannelNums,2) = i
              sawWhat = sawChannelColon
            else if ( next == sawChannelplus ) then
              if ( i <= myChannelNumbers(numChannelNums,2) ) &
                & call announce_error ( channelsDecrease, signal_string, &
                  & tree_index )
              numChannelNums = numChannelNums + 1
              myChannelNumbers(numChannelNums,:) = i
              sawWhat = sawChannel
            end if
          end do
        else
          call announce_error ( notLabel, signal_string, tree_index, &
            & expected = spec_indices((/ s_band, s_radiometer, &
            &                            s_spectrometerType /)) )
        end if
      else
        if ( verify(tokenText(1:1), 'Bb') == 0 ) then
          ! Maybe we need to splice out the sideband identifier.
          ! First, try it as-is.
          decl = get_decl ( the_token%string_index, type=label )
          if ( decl%type /= label ) then
            mySideband = scan(tokenText,'UuLl')
            if ( mySideband /= 0 ) then
              the_token%string_index = &
                & enter_terminal(tokenText(:mySideband-1) // &
                & trim(tokenText(mySideband+1:)), t_identifier)
              if ( verify(tokenText(mySideband:mySideband),'Uu') == 0 ) then
                mySideband = 1
              else
                mySideband = -1
              end if
            end if
          end if
        end if
        decl = get_decl ( the_token%string_index, type=label )
        if ( decl%type /= label ) then
          call announce_error ( notLabel, signal_string, tree_index, &
            & expected = spec_indices((/ s_band, s_radiometer, &
            &                            s_spectrometerType /)) )
        else
          select case ( get_spec_id(decl%tree) )
          case ( s_band )
            next = sawBand
            band_i = decoration(decl%tree)
          case ( s_radiometer )
            next = sawRadiometer
            radiometer_i = decoration(decl%tree)
          case ( s_spectrometerType )
            next = sawSpectrometerType
            spectrometer_i = decoration(decl%tree)
          case default
            call announce_error ( notLabel, signal_string, tree_index, &
              & expected = spec_indices((/ s_band, s_radiometer, &
              &                            s_spectrometerType /)) )
          end select
        end if
      end if

      ! Verify they've come in the correct order
      if ( error == 0 .and. next <= sawWhat ) &
        & call announce_error ( outOfOrder, signal_string, tree_index )
      sawWhat = next

      ! Check for a suffix; if there is one, verify that it's OK.
      call lex_signal ( signal_string, where, the_token )
      if ( the_token%class == t_colon ) then
        call lex_signal ( signal_string, where, the_token )
        if ( the_token%class /= t_identifier ) then
          call announce_error ( unexpected, signal_string, tree_index, &
            & expected = (/ t_identifier /) )
    exit
        end if
        call get_string ( the_token%string_index, tokenText, .true., .true. )
        select case ( next )
        case ( sawBand )           ! Band:Suffix
          call get_string ( bands(band_i)%suffix, testText, .true., .true. )
          if ( testText /= tokenText ) &
            & call announce_error ( wrongSuffix, signal_string, tree_index, &
              & expected = (/ the_token%string_index /) )
        case ( sawRadiometer )     ! Radiometer:Suffix
          call get_string ( radiometers(radiometer_i)%suffix, testText, &
            & .true., .true. )
          if ( testText /= tokenText ) &
            & call announce_error ( wrongSuffix, signal_string, tree_index, &
              & expected = (/ the_token%string_index /) )
        case default
          call announce_error ( noSuffix, signal_string, tree_index )
        end select
        call lex_signal ( signal_string, where, the_token )
      else if ( the_token%class == t_minus ) then
        call lex_signal ( signal_string, where, the_token )
        if ( the_token%class /= t_identifier ) then
          call announce_error ( unexpected, signal_string, tree_index, &
            & expected = (/ t_identifier /) )
        else
          if ( next /= sawSpectrometerType ) &
            & call announce_error ( noSuffix, signal_string, tree_index )
          call get_string ( the_token%string_index, tokenText, strip=.true. )
          read ( tokenText, *, iostat=status ) spectrometer_n
          if ( status /= 0 ) &
            & call announce_error ( notNumber, signal_string, tree_index )
        end if
        call lex_signal ( signal_string, where, the_token )
      end if
      if ( the_token%class == t_end_of_input ) &
    exit
      if ( the_token%class /= t_dot ) &
        & call announce_error ( unexpected, signal_string, tree_index, &
          & expected = (/ t_dot /) )
    end do o

    if ( error == 0 .and. numChannelNums > 0 ) then
      call allocate_test ( myChannels, &
        & maxval(myChannelNumbers(:numChannelNums,:)), &
        & "MyChannels", moduleName, &
        & lowBound=minval(myChannelNumbers(:numChannelNums,:)) )
      myChannels = .false.
      do i = 1, numChannelNums
        myChannels(myChannelNumbers(i,1):myChannelNumbers(i,2)) = .true.
      end do
    end if
    ! Now we have all of the pieces assembled and identified.  Find the
    ! set of signals that satisfies the specified criteria.
    if ( radiometer_i < 0 .and. band_i < 0 .and. spectrometer_i < 0 ) &
      call announce_error ( notEnough, signal_string, tree_index )
    if ( error == 0 ) call scanSignals

    if ( present(sideband) ) sideband = mySideband
    if ( present(channels) ) then
      call deallocate_test ( channels, "Channels", moduleName )
      channels => myChannels
      nullify ( myChannels )
    end if
    call deallocate_test ( myChannels, "MyChannels", moduleName )

  contains

    ! -------------------------------------------  Announce_Error  -----
    subroutine Announce_Error ( code, string, tree, expected )
      integer, intent(in) :: Code         ! Code index for the error
      character(len=*), intent(in) :: String   ! The string being analyzed
      integer, intent(in), optional :: Tree         ! Tree index
      integer, intent(in), optional :: Expected(:)  ! Expected tokens

      integer :: I
      integer :: myTree

      error = max(error,1)
      myTree = 0
      if ( present(tree) ) myTree = tree
      call output ( '***** At ' )
      if ( mytree > 0 ) then
        call print_source ( source_ref(mytree) )
      else
        call output ( '(no lcf tree available)' )
      end if
      call output ( ', ' )
      call output ( 'In column ' )
      call output ( where )
      call output ( ' of ' )
      call output ( trim(string) )
      call output ( ', ' )
      select case ( code )
      case ( channelsDecrease )
        call output ( 'the channel numbers are not monotone nondecreasing', &
          & advance='yes' )
      case ( invalid )
        call output ( 'does not specify a valid signal.', advance='yes' )
      case ( malformed )
        call output ( 'expected digits after the first letter.', advance='yes' )
      case ( outOfOrder )
        call output ( 'the token is out-of-order.', advance='yes' )
      case ( noSuffix )
        call output ( 'a suffix is not permitted at this point.', advance='yes' )
      case ( notEnough )
        call output ( 'not enough information given.', advance='yes' )
      case ( notLabel )
        call output ( 'the token is not the label of a ' )
        do i = 1, size(expected)
          call output ( '"' )
          call display_string ( expected(i) )
          call output ( '", ' )
        end do
        call output ( '"switch" or "channel".', advance='yes' )
      case ( notNumber )
        call output ( 'expected a number.', advance='yes' )
      case ( unexpected )
        call output ( 'the input was incorrect.  Expected ' )
        do i = 1, size(expected)
          call dump_symbol_class ( expected(i) )
          if ( size(expected) > 1 ) then
            if ( i < size(expected) - 1 ) call output ( ', ' )
            if ( i == size(expected) - 1 ) call output ( ' or ' )
          end if
        end do
        call output ( '.', advance='yes' )
      case ( wrongSuffix )
        call output ('the wrong suffix is specified.', advance='yes' )
      case default
        call output ( 'No code for error code ' )
        call output ( code )
        call output ( ' in Parse_signal_m%Announce_Error', advance='yes' )
      end select
    end subroutine Announce_Error

    ! ----------------------------------------------  ScanSignals  -----
    subroutine ScanSignals
      ! Scan the signals database for entries that satisfy the criteria
      ! given by the tokens in the signal string.
      ! Allocate Signal_Indices and fill the array.

      logical :: BandMatch(size(signals))
      logical :: ChannelMatch(size(signals))
      integer :: HowMany
      integer :: I
      logical :: RadiometerMatch(size(signals))
      integer :: Spectrometer           ! Index from the signals database
      logical :: SpectrometerMatch(size(signals))
      logical :: SpectrometerTypeMatch(size(signals))
      logical :: SwitchMatch(size(signals))

      ! Initialize
      bandMatch = band_i < 0
      channelMatch = numChannelNums <= 0
      radiometerMatch = radiometer_i < 0
      spectrometerMatch = spectrometer_n < 0
      spectrometerTypeMatch = spectrometer_i < 0
      switchMatch = switch < 0

      do i = 1, size(signals)
        spectrometer = signals(i)%spectrometerType
        if ( band_i >= 0 ) bandMatch(i) = band_i == signals(i)%band
        if ( numChannelNums > 0 ) &
          & channelMatch(i) = &
            & lbound(myChannels,1) >= lbound(signals(i)%frequencies,1) .and. &
            & ubound(myChannels,1) <= ubound(signals(i)%frequencies,1)
        if ( radiometer_i >= 0 ) radiometerMatch(i) = &
          & radiometer_i == signals(i)%radiometer
        if ( spectrometer_n >= 0 ) spectrometerMatch(i) = &
          & spectrometer_n ==signals(i)%spectrometer
        if ( spectrometer_i >= 0 ) spectrometerTypeMatch(i) = &
          & spectrometer_i == spectrometer
        if ( switch >= 0 ) switchMatch(i) = switch ==signals(i)%switch
      end do
      howMany = count ( bandMatch .and. channelMatch .and. radiometerMatch &
        & .and. spectrometerMatch .and. SpectrometerTypeMatch .and. &
        & switchMatch )
      if ( present(onlyCountEm) ) onlyCountEm = howMany
      if ( howMany == 0 ) then
        where = 1
        call announce_error ( invalid, signal_string, tree_index )
      else if ( .not. present(onlyCountEm) ) then
        call allocate_test ( signal_indices, howMany, moduleName, &
          & 'signal_indices' )
        howMany = 0
        do i = 1, size(signals)
          if ( bandMatch(i) .and. channelMatch(i) .and. radiometerMatch(i) &
          & .and. spectrometerMatch(i) .and. SpectrometerTypeMatch(i) .and. &
          & switchMatch(i) ) then
            howMany = howMany + 1
            signal_indices(howMany) = i
          end if
        end do
      end if
    end subroutine ScanSignals

  end subroutine Parse_Signal

  logical function not_used_here()
    !---------------------------- RCS Ident Info -------------------------------
    character (len=*), parameter :: IdParm = &
      & "$Id$"
    character (len=len(idParm)) :: Id = idParm
    !---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Parse_Signal_M

! $Log$
! Revision 2.22  2005/05/02 22:58:24  vsnyder
! Modify interface for Get_Individual_Signals
!
! Revision 2.21  2005/03/24 01:38:36  vsnyder
! Spiff up an error message
!
! Revision 2.20  2005/02/05 00:05:07  vsnyder
! Correct the low bound for channels in Expand_signal_list.  Handle
! disassociated channels argument from parse_signal correctly in
! Get_Individual_Signals.  Move CVS Id into not_used_here.  Some
! cannonball polishing.
!
! Revision 2.19  2005/01/12 23:58:47  vsnyder
! All in Get_Individual_Sgnals: Detect DSB signals correctly.  Correct an
! indexing error.  Include the channel.
!
! Revision 2.18  2005/01/12 03:08:38  vsnyder
! Add Get_Individual_Signals
!
! Revision 2.17  2004/03/24 22:45:59  livesey
! Increased buffer sizes
!
! Revision 2.16  2003/05/05 23:00:05  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.15  2003/03/07 03:17:28  livesey
! Added Expand_Signal_List
!
! Revision 2.14.2.1  2003/03/01 02:40:04  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.14  2003/01/25 04:14:13  vsnyder
! Get rid of USEs for stuff not actually used
!
! Revision 2.13  2002/11/06 00:13:32  pwagner
! Should not dump core if parse_signal called with tree_index 0
!
! Revision 2.12  2002/10/08 00:09:13  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.11  2001/06/07 21:59:41  pwagner
! Added Copyright statement
!
! Revision 2.10  2001/04/26 02:33:03  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.9  2001/04/13 20:58:48  vsnyder
! Add 'OnlyCountEm' argument
!
! Revision 2.8  2001/04/11 20:19:46  vsnyder
! Look for channels in signals instead of spectrometers
!
! Revision 2.7  2001/04/11 02:10:45  vsnyder
! Fix a typo
!
! Revision 2.6  2001/04/10 17:59:54  vsnyder
! Remove sideband field from signal
!
! Revision 2.5  2001/04/09 20:15:53  vsnyder
! Tighter bound on myChannelNumbers
!
! Revision 2.4  2001/04/06 20:15:36  vsnyder
! Implement syntax for specifying several channels
!
! Revision 2.3  2001/03/16 00:34:50  vsnyder
! Correct handling of the Band database
!
! Revision 2.2  2001/03/15 23:57:27  vsnyder
! OOPS, forgot to put 'Log' at the end
!
