! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Dump_Options

  ! Options for Dump_0, Dump_1, Diff_1, and maybe others.
  ! Variables shared between those modules and derived from the options.


! Most of the optional parameters have default values
! logically set to FALSE or 0, ' ',  or '*' where appropriate

!  optional args
!  (dumps and diffs if the same)
!      arg            meaning                                     default
!      ---            -------                                     -------
!    fillvalue        skip dumping lines containg only fillValues    0
!    width            how many values printed per line            depends
!    format           fortran format used to print                depends
!    lbound           lower bound of 1st index                       1
!    options          (see below)                                    ''

!  (diffs)
!      arg            meaning                                     default
!      ---            -------                                     -------
!    fillvalue        don't diff where array elements are this    -999.99

! The format optional arg defaults to SDFORMATDEFAULT for floating pt. arrys
! For integer arrays it defaults to i6 or i_INTPLACES_
! For complex arrays it defaults to SDFORMATDEFAULTCMPLX
! If set to '(*)', for floating point and complex arrays it
! will be the least number of spaces wide enough to contain the
! largest array element printed according to the default format
! If set to '(*.m)', for floating point and complex arrays it
! will be the least number of spaces wide enough to contain the
! largest array element with m spaces after the decimal point

! The meaning of options has replaced the older logical arguments
! if the options is present and contains the following characters:
! (for dump or diff)
!   character         meaning
!      ---            -------
!       @              show differences in the goldbrick-style
!       B              show Bandwidth, % of array that is non-zero
!       H              show rank, shape of array
!       L              laconic; skip printing name, size of array
!       N              Show where NaNs and Infs are located
!       R              rms       -- min, max, etc.
!       b              table of % vs. amount of differences (pdf)
!       c              clean
!       g              gaps      
!       l              collapse (last index)
!       r              ratios    -- min, max, etc. of differences' ratios
!       s              stats     -- number, % of differences
!       p              transpose 
!       t              trim      
!       u              unique    
!       v              verbose
!       w              wholearray
!       W[i]           wholearray, looping over ith index (for rank 3 and 4 arrays only)
!       1 or 2 or ..   ignored; calling routine is free to interpret

! An exception is the behavior of wholearray:
! if all {HRblrs} are FALSE, i.e. unset, the whole array is dumped (or diffed)
! if any is TRUE the whole array will be dumped only if
! w or wholearray is set to TRUE

! in the above, a string list is a string of elements (usu. comma-separated)
!
! For an example of the goldbrick style of showing differences,
!  ** Reference values (min : max, rms): -343.791 : 4.67836, 3.52708
!  ** Max. absolute: 0.000439644 ( = 4.39644e+16 fractional )
!  ** Max. fractional: 4.39644e+16 ( = 0.000439644 absolute )


  implicit none
  public

  save

  character, public, parameter :: AfterSub = '#'
  character(len=*), parameter :: DefaultPCTFormat = '(0pf6.1)'

  ! These are the possible option characters to dumps, diffs
  character, parameter :: Dopt_AuBrick     = '@'
  character, parameter :: Dopt_Bandwidth   = 'B'
  character, parameter :: Dopt_Clean       = 'c'
  character, parameter :: Dopt_Collapse    = 'l'
  character, parameter :: Dopt_Cyclic      = 'y'
  character, parameter :: Dopt_Direct      = 'd'
  character, parameter :: Dopt_Gaps        = 'g'
  character, parameter :: Dopt_Laconic     = 'L'
  character, parameter :: Dopt_NaNs        = 'N'
  character, parameter :: Dopt_OnlyWholeArray = '' ! Set from other options
  character, parameter :: Dopt_Ratios      = 'r'
  character, parameter :: Dopt_Rms         = 'R'
  character, parameter :: Dopt_Verbose     = 'v'
  character, parameter :: Dopt_Shape       = 'H'
  character, parameter :: Dopt_Stats       = 's'
  character, parameter :: Dopt_Table       = 'b'
  character, parameter :: Dopt_Transpose   = 'p'
  character, parameter :: Dopt_Trim        = 't'
  character, parameter :: Dopt_Unique      = 'u'
  character, parameter :: Dopt_WholeArray  = 'w'

  ! These are the indices of those options in the Dopts array
  enum, bind(c)
    enumerator :: AuBrick = 1
    enumerator :: Bandwidth
    enumerator :: Clean
    enumerator :: CollapseIt
    enumerator :: Cyclic
    enumerator :: Direct
    enumerator :: ItsShape
    enumerator :: Gaps
    enumerator :: Laconic
    enumerator :: NaNs
    enumerator :: OnlyWholeArray
    enumerator :: Ratios
    enumerator :: Rms
    enumerator :: Verbose
    enumerator :: Stats
    enumerator :: Table
    enumerator :: Transpose
    enumerator :: TrimIt
    enumerator :: Unique
    enumerator :: WholeArray
  end enum

  integer, private, parameter :: NumOpt = WholeArray

  logical :: NameHasBeenPrinted = .false.

!     (parameters)
! CollapseOptions          options determining what and how to dump collapsed
!                           representations of multidimensional arrays
! DefaultDiffOptions       switches to set default DIFF values for CLEAN, TRIM, etc.
! DefaultDumpOptions       same as above, but for DUMP
! DiffRMSMeansRMS          print abs min, max, etc. when DIFF has RMS set TRUE
! DontDumpIfAllEqual       don't dump every element of a constant array
! DumpTableSide            what side to place headers when dumping tables
! FilterFillsFromRMS       exclude fill values when calculating rms, etc.
!                           (not implemented yet)
! IntPlaces                How many places to print when dumping integer values
! MaxNumNANS               How many NaNs can we show where they are
! NameOnEachLine           item name to print on each output line
! PCTFormat                use this format to print % with '-s' diff option
! RMSFormat                use this format to print min, max, rms, etc.
! SDFormatDefault          use this format to print s.p., d.p. by default
! StatsOnOneLine           stats, rms each printed on a single line

  ! The following character strings can include one or more options listed above
  ! E.g., '-crt' turns on Clean, RMS, and TrimIt
  character(len=8)  :: DefaultDiffOptions = ' '
  character(len=8)  :: DefaultDumpOptions = ' '

  integer           :: DefaultMaxLon      = 128
  integer           :: DefaultWidth       = 10
  character(len=8)  :: DumpTableSide      = 'top'
  logical           :: DiffRMSMeansRMS    = .false.
  logical           :: DontDumpIfAllEqual = .true.
  logical           :: FilterFillsFromRMS = .false.
  logical           :: PrintFillValue     = .true.
  logical           :: PrintNameIfDiff    = .true.
  logical           :: StatsOnOneLine     = .true.
  character(len=16) :: NameOnEachLine     = ' '

  ! This determines how a higher-rank array is collapsed to a lower-rank one
  character(len=16) :: CollapseOptions = 'num[+]all[+]'

  ! These determine how dumped numerical data (s.p. or d.p.) will be formatted
  character(len=2)  :: IntPlaces = '6' ! how many places
  integer           :: MaxNumNANs= 60  ! how many NaNs to show
  character(len=16) :: PCTFormat = '*' ! * means default format
  character(len=16) :: RMSFormat = '*' ! * means default format
  character(len=16) :: SDFormatDefault = '(1pg14.6)'
  character(*), parameter :: SDFormatDefaultCmplx = &
    & '(1x,"(",1pg13.6,",",1pg13.6,")")'

  type :: Option_T
    character(10) :: Name  ! Option's name
    character :: Char      ! Letter that selects option
    logical :: V           ! Value of the option
  end type Option_T

  ! Options for Dumps and Diffs
  type(option_t) :: Dopts(numOpt)

  data dopts(aubrick)    / option_t('aubrick',    dopt_aubrick,    .false. ) /
  data dopts(bandwidth)  / option_t('bandwidth',  dopt_bandwidth,  .false. ) /
  data dopts(clean)      / option_t('clean',      dopt_clean,      .false. ) /
  data dopts(collapseIt) / option_t('collapse',   dopt_collapse,   .false. ) /
  data dopts(cyclic)     / option_t('cyclic',     dopt_cyclic,     .false. ) /
  data dopts(direct)     / option_t('direct',     dopt_direct,     .false. ) /
  data dopts(gaps)       / option_t('gaps',       dopt_gaps,       .false. ) /
  data dopts(itsShape)   / option_t('shape',      dopt_shape,      .false. ) /
  data dopts(laconic)    / option_t('laconic',    dopt_laconic,    .false. ) /
  data dopts(NaNs)       / option_t('NaNs',       dopt_nans,       .false. ) /
  data dopts(onlyWholeArray) / option_t('OnlyWholeArray', '',      .false. ) /
  data dopts(ratios)     / option_t('ratios',     dopt_ratios,     .false. ) /
  data dopts(RMS)        / option_t('RMS',        dopt_rms,        .false. ) /
  data dopts(stats)      / option_t('stats',      dopt_stats,      .false. ) /
  data dopts(table)      / option_t('table',      dopt_transpose,  .false. ) /
  data dopts(transpose)  / option_t('transpose',  dopt_transpose,  .false. ) /
  data dopts(trimIt)     / option_t('trim',       dopt_trim,       .false. ) /
  data dopts(unique)     / option_t('unique',     dopt_unique,     .false. ) /
  data dopts(verbose)    / option_t('verbose',    dopt_verbose,    .false. ) /
  data dopts(wholeArray) / option_t('wholeArray', dopt_wholeArray, .false. ) /

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains
  
  ! ---------------------------------------------- DumpDumpOptions -----
  subroutine DumpDumpOptions( OptionString, ThisIsADiff )
    ! Show 
    ! (if current OptionString supplied) actual dump, diff options
    ! (if arg is "?") available options and their meanings
    ! (if no arg) default options
    use HighOutput, only: OutputNamedValue
    use MLSStrings, only: Trim_safe
    use Output_m, only: Blanks, Output
    character(len=*), intent(in), optional :: OptionString
    logical, intent(in), optional          :: ThisIsADiff
    character(len=1), parameter :: FillChar = '1' ! fill blanks with '. .'
    integer :: I
    logical :: MyDiff
    myDiff = .false.
    if ( present(thisIsADiff) ) myDiff = thisIsADiff
    if ( .not. present(OptionString) ) then
      call blanks(90, fillChar='-', advance='yes')
      call output(' ------------------------ Summary of automatic Dump, Diff options'      , advance='no')
      call output(' ------------------------ ', advance='yes')
      call outputNamedValue ( 'character printed between row, col id and data', aftersub, advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      call outputNamedValue ( 'default DIFF switches for CLEAN, TRIM, etc.', trim_safe(DefaultDiffOptions), advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      call outputNamedValue ( 'default DUMP switches for CLEAN, TRIM, etc.', trim_safe(DefaultDumpOptions), advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      call outputNamedValue ( 'print abs min, max, etc. when DIFF has RMS set TRUE?', DiffRMSMeansRMS, advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      call outputNamedValue ( 'skip dumping every element of a constant array?', DiffRMSMeansRMS, advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      call outputNamedValue ( 'what side to place headers when dumping tables', trim(DumpTableSide), advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      call outputNamedValue ( 'print stats all on one line?', StatsOnOneLine, advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      call outputNamedValue ( 'pct output format', trim_safe(PCTFormat), advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      call outputNamedValue ( 'rms output format', trim_safe(RMSFormat), advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      call outputNamedValue ( 'numeric output format', trim_safe(SDFormatDefault), advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      call outputNamedValue ( 'complex output format', trim_safe(sdFormatDefaultCmplx), advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      call blanks(90, fillChar='-', advance='yes')
      do i = 1, size(dopts)
        call outputNamedValue ( dopts(i)%name//'?', dopts(i)%v, advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      end do
    else if ( index( optionString, "?" ) > 0 ) then
      call output( 'The meaning of options is determined by the presence or absence', advance='yes' )
      call output( 'of specific characters in the options string', advance='yes' )
      call output( 'If options is present and contains the following characters:', advance='yes' )
      call output( '(for dump or diff)', advance='yes' )
      call output( '  character         meaning', advance='yes' )
      call output( '     ---            -------', advance='yes' )
      ! The following are in order according to their option letters
      call output( '     ' // dopt_aubrick     // '              show aubrick, diffs in goldbrick style', advance='yes' )
      call output( '     ' // dopt_bandwidth   // '              show Bandwidth, % of array that is non-zero', advance='yes' )
      call output( '     ' // dopt_shape       // '              show rank, shape of array', advance='yes' )
      call output( '     ' // dopt_laconic     // '              laconic; skip printing name, size of array', advance='yes' )
      call output( '     ' // dopt_NaNs        // '              show where NaNs and Infs are located', advance='yes' )
      call output( '     ' // dopt_RMS         // '              rms       -- min, max, etc.', advance='yes' )
      call output( '     ' // dopt_table       // '              table of % vs. amount of differences (pdf)', advance='yes' )
      call output( '     ' // dopt_clean       // '              clean', advance='yes' )
      call output( '     ' // dopt_direct      // '              direct', advance='yes' )
      call output( '     ' // dopt_gaps        // '              gaps      ', advance='yes' )
      call output( '     ' // dopt_collapse    // '              collapse (last index)', advance='yes' )
      call output( '     ' // dopt_ratios      // '              ratios    -- min, max, etc. of difference ratios', advance='yes' )
      call output( '     ' // dopt_stats       // '              stats     -- number, % of differences', advance='yes' )
      call output( '     ' // dopt_transpose   // '              transpose (for rank 2 arrays only)', advance='yes' )
      call output( '     ' // dopt_trim        // '              trim      ', advance='yes' )
      call output( '     ' // dopt_unique      // '              unique    ', advance='yes' )
      call output( '     ' // dopt_verbose     // '              verbose   ', advance='yes' )
      call output( '     ' // dopt_cyclic      // '              cyclic    ', advance='yes' )
      call output( '     ' // dopt_wholeArray  // '              wholearray', advance='yes' )
      call output( '     W[i]           wholearray, looping over ith index', advance='yes' )
      call output( '                    (for rank 3 and 4 arrays only)', advance='yes' )
      call output( '     1 or 2 or ..   ignored; calling routine is free to interpret', advance='yes' )
      call output( ' ', advance='yes' )
      call output( 'An exception is the behavior of w (wholearray):', advance='yes' )
      call output( 'if all {@HNRblrs} are FALSE, i.e. unset, the whole array is dumped (or diffed)', advance='yes' )
      call output( 'if any are TRUE the whole array will be dumped only if', advance='yes' )
      call output( 'w is present or wholearray is set to TRUE', advance='yes' )
    else
      call outputNamedValue ( 'options', trim_safe(optionString), advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      call outputNamedValue ( 'thisIsADiff?', myDiff, advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      do i = 1, size(dopts)
        call outputNamedValue ( dopts(i)%name//'?', dopts(i)%v, advance='yes', &
        & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
      end do
    end if
  end subroutine DumpDumpOptions

  ! ------------------------------------------  RestoreDumpConfig  -----
  ! Restore default values for dump settings
  subroutine RestoreDumpConfig
    DefaultDiffOptions        = ' '
    DefaultDumpOptions        = ' '
    DefaultMaxLon             = 128
    DefaultWidth              = 10
    DumpTableSide             = 'top'
    DiffRMSMeansRMS           = .false.
    DontDumpIfAllEqual        = .true.
    FilterFillsFromRMS        = .false.
    printFillValue            = .true.
    printNameIfDiff           = .true.
    StatsOnOneLine            = .true.
    collapseOptions           = 'num[+]all[+]'
    IntPlaces                 = '6' ! how many places
    PCTFormat                 = '*' ! * means default format
    RMSFormat                 = '*' ! * means default format
    SDFormatDefault           = '(1pg14.6)'
  end subroutine RestoreDumpConfig

  ! ------------------------------------------------  Set_Options  -----
  subroutine Set_Options ( OptionString, DefaultOptionString )
    ! Set the values in Options according to what characters are in
    ! DefaultOptionString, if it's present, else set them false.  Then
    ! turn on options in OptionString.
    character(*), intent(in), optional :: OptionString
    character(*), intent(in), optional :: DefaultOptionString
    integer :: I

    dopts%v = .false.
    if ( present(defaultOptionString) ) then
      do i = 1, size(dopts)
        dopts(i)%v = index(defaultOptionString, dopts(i)%char) > 0
      end do
    else
      dopts%v = .false.
    end if
    if ( present(optionString) ) then
      do i = 1, size(dopts)
        dopts(i)%v = dopts(i)%v .or. index(optionString, dopts(i)%char) > 0
      end do
    end if

  end subroutine Set_Options

  ! -----------------------------------------------  TheDumpBegins -----
  ! A warm-up subroutine that transfers the options string to the Dopts
  subroutine TheDumpBegins ( Options, ThisIsADiff )
    use Output_m, only: StampOptions
    character(len=*), intent(in), optional :: Options
    logical, intent(in), optional          :: ThisIsADiff
    logical :: MyDiff
    ! Executable
    ! Were we called with the trigger '?'?
    if ( present(options) ) then
      if ( index( options, '?' ) > 0 ) then
        call DumpDumpOptions ( options )
        stop
      endif
    endif
    myDiff = .false.
    if ( present(thisIsADiff) ) myDiff = thisIsADiff
    nameHasBeenPrinted = .false.
    stampOptions%neverStamp = .true. ! So we don't interrupt tables of numbers
    call set_options ( options, &
      & merge(defaultDiffOptions,defaultDumpOptions,myDiff) )
    dopts(wholeArray)%v = dopts(wholeArray)%v .or. &
      & .not. ( dopts(AuBrick)%v .or. dopts(bandwidth)%v .or. dopts(collapseIt)%v .or. &
      &         dopts(ratios)%v .or. dopts(RMS)%v .or. dopts(itsShape)%v .or. &
      &         dopts(stats)%v .or. dopts(table)%v .or. dopts(NaNs)%v )
    dopts(onlyWholeArray)%v = dopts(wholeArray)%v .and. &
      & .not. ( dopts(AuBrick)%v .or. dopts(bandwidth)%v .or. dopts(collapseIt)%v .or. &
      &         dopts(ratios)%v .or. dopts(RMS)%v .or. dopts(itsShape)%v .or. &
      &         dopts(stats)%v .or. dopts(table)%v .or. dopts(NaNs)%v)
    nameHasBeenPrinted = nameHasBeenPrinted .or. dopts(laconic)%v
  end subroutine TheDumpBegins

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Dump_Options

! $Log$
! Revision 2.4  2016/12/14 01:20:55  pwagner
! Print meaning of all options if options is '?'
!
! Revision 2.3  2016/09/09 20:09:51  pwagner
! Added Au (Gold) brick option removing some hay from the stack of statistics
!
! Revision 2.2  2016/07/28 03:29:28  vsnyder
! Moved a bunch of comments here from Dump_0.  Repaired typo that confused
! "Clean" option with "Collape" option.
!
! Revision 2.1  2016/07/28 01:41:48  vsnyder
! Initial Commit
!
