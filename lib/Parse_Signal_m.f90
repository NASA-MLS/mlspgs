module Parse_Signal_M

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Declaration_Table, only: Decls, Get_Decl, Label
  use Init_MLSSignals_m, only: S_Band, S_Module, S_Radiometer, S_Signal, &
    & S_SpectrometerType, Spec_First
  use Lexer_Core, only: Print_Source, Token
  use Lexer_m, only: Lex_Signal
  use MLSSignals_m, only: Bands, Radiometers, Signals, SpectrometerTypes
  use MoreTree, only: Get_Spec_ID
  use Output_m, only: Output
  use String_table, only: Display_string, Get_String
  use Symbol_Table, only: Dump_symbol_class
  use Symbol_Types, only: T_Colon, T_Dot, T_End_of_input, T_Identifier, &
    & T_Minus
  use Tree, only: Decoration, Source_Ref

  implicit NONE
  private
  public :: Parse_Signal

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  subroutine Parse_Signal ( Signal_String, Signal_Indices, Spec_indices, &
    & Tree_Index )

  ! Parse a signal string.  Return the index of the signal in the signal
  ! database.

    character(len=*), intent(in) :: Signal_String ! Input
    integer, pointer :: Signal_Indices(:)         ! In the signals database
    ! They are allocated here using Allocate_Test, so don't start with an
    ! undefined pointer!  Upon return, if Signal_Indices is not associated,
    ! an error occurred.
    integer, intent(in) :: Spec_indices(spec_first:)   ! Indices in string
    !                                               table; For error messages
    integer, intent(in), optional :: Tree_Index   ! To get line and column
    !                                               numbers for messages

    integer :: Band_i         ! Database index
    integer :: Channel        ! Channel number, from the token text
    type(decls) :: Decl       ! Declaration of an identifier
    integer :: Error          ! >0 indicates an error
    integer :: I              ! Temporary, loop inductor, subscript
    integer :: Next           ! One of Saw... parameters below, to indicate
    !                           the category of the next identifier
    integer :: Radiometer_i   ! Database index
    integer :: SawWhat        ! What has been seen so far.  See Saw... below
    integer :: Spectrometer_i ! Database index
    integer :: Spectrometer_n ! Number after "-", in e.g. FB25-3
    integer :: Status         ! From an internal READ
    integer :: Switch         ! Switch number, from the token text
    character(len=5) :: TestText   ! From a "%suffix" field via sub_rosa
    type(token) :: The_Token  ! that results from lexing the signal string
    character(len=5) :: TokenText  ! Duh
    integer :: Where     ! The current analysis point in the signal string

    ! States.  The pieces of the radiometer specification must be presented
    ! in order, but pieces can be missing.
    integer, parameter :: SawNothing = 0
    integer, parameter :: SawRadiometer = sawNothing + 1
    integer, parameter :: SawBand = sawRadiometer + 1
    integer, parameter :: SawSwitch = sawBand + 1
    integer, parameter :: SawSpectrometerType = sawSwitch + 1
    integer, parameter :: SawChannel = sawSpectrometerType + 1

    ! Parameters for error message codes
    integer, parameter :: Invalid = 1                ! Invalid combination
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
    channel = -1
    error = 0
    next = sawNothing
    radiometer_i = -1
    sawWhat = sawNothing
    call deallocate_test ( signal_indices, moduleName, "Signal_Indices" )
    spectrometer_i = -1
    spectrometer_n = -1
    switch = -1
    where = 0

    ! Get the tokens; verify they are the correct kinds of thing-o's.
    do
      call lex_signal ( signal_string, where, the_token )
      if ( the_token%class == t_end_of_input ) &
    exit
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
          channel = i
        else
          call announce_error ( notLabel, signal_string, tree_index, &
            & expected = spec_indices((/ s_band, s_radiometer, &
            &                            s_spectrometerType /)) )
        end if
      else
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

      ! Check for a suffix; if there is one, verify that it's OK.
      call lex_signal ( signal_string, where, the_token )
      if ( the_token%class == t_colon ) then
        call lex_signal ( signal_string, where, the_token )
        if ( the_token%class /= t_identifier ) &
          & call announce_error ( unexpected, signal_string, tree_index, &
            & expected = (/ t_identifier /) )
        call get_string ( the_token%string_index, tokenText, .true., .true. )
        select case ( next )
        case ( sawBand )
          call get_string ( bands(band_i)%suffix, testText, .true., .true. )
          if ( testText /= tokenText ) &
            & call announce_error ( wrongSuffix, signal_string, tree_index, &
              & expected = (/ the_token%string_index /) )
        case ( sawRadiometer )
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
    end do

    ! Now we have all of the pieces assembled and identified.  Find the
    ! set of signals that satisfy the specified criteria.
    if ( radiometer_i < 0 .and. band_i < 0 .and. spectrometer_i < 0 ) &
      call announce_error ( notEnough, signal_string, tree_index )
    if ( error == 0 ) call scanSignals

  contains

    ! -------------------------------------------  Announce_Error  -----
    subroutine Announce_Error ( code, string, tree, expected )
      integer, intent(in) :: Code         ! Code index for the error
      character(len=*), intent(in) :: String   ! The string being analyzed
      integer, intent(in), optional :: Tree         ! Tree index
      integer, intent(in), optional :: Expected(:)  ! Expected tokens

      integer :: I

      error = max(error,1)
      if ( present(tree) ) then
        call output ( 'At ' )
        call print_source ( source_ref(tree) )
        call output ( ', ' )
      end if
      call output ( 'In column ' )
      call output ( where )
      call output ( ' of ' )
      call output ( trim(string) )
      call output ( ', ' )
      select case ( code )
      case ( invalid )
        call output ( 'does not specify a valid signal.', 'yes' )
      case ( malformed )
        call output ( 'expected digits after the first letter.', 'yes' )
      case ( outOfOrder )
        call output ( 'the token is out-of-order.', 'yes' )
      case ( noSuffix )
        call output ( 'a suffix is not permitted at this point.', 'yes' )
      case ( notEnough )
        call output ( 'not enough information given.', 'yes' )
      case ( notLabel )
        call output ( 'the token is not the label of a ' )
        do i = 1, size(expected)
          call output ( '"' )
          call display_string ( expected(i) )
          call output ( '", ' )
        end do
        call output ( '"switch" or "channel".', advance = 'yes' )
      case ( notNumber )
        call output ( 'expected a number.', 'yes' )
      case ( unexpected )
        call output ( 'the input was incorrect.  Expected ' )
        do i = 1, size(expected)
          call dump_symbol_class ( expected(i) )
          if ( size(expected) > 1 ) then
            if ( i < size(expected) - 1 ) call output ( ', ' )
            if ( i == size(expected) - 1 ) call output ( ' or ' )
          end if
        end do
        call output ( '.', advance = 'yes' )
      end select
    end subroutine Announce_Error

    ! ----------------------------------------------  ScanSignals  -----
    subroutine ScanSignals
      ! Scan the signals database for entries that satisfy the criteria
      ! given by the tokens in the signal string.
      ! Allocate Signal_Indices , and fill the array.

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
      channelMatch = channel < 0
      radiometerMatch = radiometer_i < 0
      spectrometerMatch = spectrometer_n < 0
      spectrometerTypeMatch = spectrometer_i < 0
      switchMatch = switch < 0

      do i = 1, size(signals)
        spectrometer = signals(i)%spectrometerType
        if ( band_i >= 0 ) bandMatch(i) = band_i == signals(i)%band
        if ( channel >= 0 ) channelMatch(i) = &
          & channel >= lbound(spectrometerTypes(spectrometer)%frequencies,1) .and. &
          & channel <= ubound(spectrometerTypes(spectrometer)%frequencies,1)
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
      if ( howMany == 0 ) then
        where = 1
        call announce_error ( invalid, signal_string, tree_index )
      else
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

end module Parse_Signal_M

! $Log$
! Revision 2.3  2001/03/16 00:34:50  vsnyder
! Correct handling of the Band database
!
! Revision 2.2  2001/03/15 23:57:27  vsnyder
! OOPS, forgot to put 'Log' at the end
!
