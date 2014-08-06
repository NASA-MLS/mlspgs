! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module DUMP_0

! Low-level dump routines -- for some arrays of intrinsic type.
! Should handle most combinations of rank and type
! Behavior depends on optional parameters
! Actual output device determined by output_m module
!
! In general, dumping a whole array of values will be presented as a matrix
! up to 10 values across
! Instead of a whole array, or in addition, one may dump a condensed summary
! showing min, max, percentages of non-zero values, etc.

! This has become too long--we may split it, putting diffs into a higher-level
! and separate dump_1.f90 module

  use BITSTUFF, only: MAXBITNUMBER, WHICHBITSARESET
  use DATES_MODULE, only: MAXUTCSTRLENGTH, &
    & REFORMATDATE, REFORMATTIME, SPLITDATETIME, TAI93S2UTC
  use HIGHOUTPUT, only: ALIGNTOFIT, BLANKSTOTAB, &
    & NUMNEEDSFORMAT, NUMTOCHARS, &
    & OUTPUTLIST, OUTPUTNAMEDVALUE, RESETTABS, SETTABS
  use IEEE_ARITHMETIC, only: IEEE_IS_FINITE
  use MLSFILLVALUES, only : BANDWIDTH, COLLAPSE, FILTERVALUES, HALFWAVES, &
    & ISFINITE, ISINFINITE, ISNAN, &
    & INFFUNCTION, NANFUNCTION, REORDERFILLVALUES, REPLACEFILLVALUES, &
    & WHEREARETHEINFS, WHEREARETHENANS
  use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_WARNING
  use MLSFINDS, only: FINDUNIQUE
  use MLSSTATS1, only: STAT_T, &
    & ALLSTATS, FILLVALUERELATION, HOWFAR, HOWNEAR, &
    & MLSMAX, MLSMEAN, MLSMIN, MLSSTDDEV, RATIOS, RESET
  use MLSSTRINGLISTS, only: CATLISTS, GETSTRINGELEMENT, NUMSTRINGELEMENTS, &
    & OPTIONDETAIL
  use MLSSTRINGS, only: DELETE, INDEXES, LOWERCASE, &
    & READINTSFROMCHARS, TRIM_SAFE, &
    & WRITEINTSTOCHARS
  use OUTPUT_M, only: OUTPUTOPTIONS, STAMPOPTIONS, &
    & BLANKS, NEWLINE, OUTPUT
  use TIME_M, only: TIME_NOW

  implicit none
  private
! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (parameters)
! AFTERSUB                 character printed between row, col id and data
! COLLAPSEOPTIONS          options determining what and how to dump collapsed
!                           representations of multidimensional arrays
! DEFAULTDIFFOPTIONS       switches to set default DIFF values for CLEAN, TRIM, etc.
! DEFAULTDUMPOPTIONS       same as above, but for DUMP
! DIFFRMSMEANSRMS          print abs min, max, etc. when DIFF has RMS set TRUE
! DONTDUMPIFALLEQUAL       don't dump every element of a constant array
! DUMPTABLESIDE            what side to place headers when dumping tables
! FILTERFILLSFROMRMS       exclude fill values when calculating rms, etc.
!                           (not implemented yet)
! INTPLACES                How many places to print when dumping integer values
! MAXNUMNANS               How many NaNs can we show where they are
! PCTFORMAT                use this format to print % with '-s' diff option
! RMSFORMAT                use this format to print min, max, rms, etc.
! SDFORMATDEFAULT          use this format to print s.p., d.p. by default
! STATSONONELINE           stats, rms each printed on a single line

!     (subroutines and functions)
! DIFF_FUN                 returns differences between scalars, arrays, etc.
! DIFF                     dump diffs between pair of arrays of numeric type
! DUMP                     dump an array to output
! DUMPDATES                dump 1-d array of tai93 (s. after 1 jan 1993)
! DUMPDUMPOPTIONS          dump module settings for dump, diff, etc.
! DUMPLISTS                dump 2-d array as a set of lists
! DUMPNAMEDVALUES          dump an array of paired names and values
! DUMPSUMS                 dump after summing successive array values
!                            ("inverse" of selfDiff)
! DUMPTABLE                dump a 2-d table of values with headers
! RESTOREDUMPCONFIG        restore default values to dump configuration
! SELFDIFF                 dump increments between successive array values
! === (end of toc) ===

! === (start of api) ===
! num diff_fun ( num value1, num value2, num auxvalue, char* options )
! diff ( array1, char* name1, array2, char* name2,
!      [fillvalue], [int width], [char* format],
!      [int lbound], [char* options] ) 
!       where array1, array2 can be 1, 2, or 3d arrays of 
!       ints, reals, or doubles, compatible in size and type
!       and fillValue is a scalar of the same type, if present
! dump ( array, char* name,
!      [fillvalue], [int width], [char* format],
!      [int lbound], [char* options] ) 
!       where array can be a 1, 2, or 3d array of
!       chars, ints, reals, or doubles,
!       and fillValue is a scalar of the same type, if present
! dump ( strlist string, char* name, [char* fillvalue], [char* options] )
! dump ( log countEmpty, strlist keys, strlist values, char* name, 
!       [char* separator], [char* options] )
! dumpDates ( dble dates(:), [int width], [char* dateFormat], [char* timeFormat] )
! dumpDumpOptions ( [char* options] )
! dumpLists ( array, char* name,
!      [int width], [char* sep],
!      [char* delims] ) 
!       where array can be a 2d array of
!       chars or ints
! dumpNamedValues ( values, strlist names,
!      [char* format, [int width], [char* options] ) 
!       where values can be a 1d array of ints or reals, and
!       names is a string list of corresponding names
! dumpSums ( array, char* name,
!      [fillvalue], [int width], [char* format],
!      [int lbound], [char* options] )
! dumpTable ( values, headers, char* headside
!      [char* format, [char* formats(:)] ) 
!       where values can be a 2d array of ints or reals, and
!       headers is an array the same size as the 2nd index of values
!       format optionally overrides the default format for the numeric type
!       formats allows you to specify a format separately column-by-column
! restoreDumpConfig
! selfdiff ( array, char* name,
!      [fillvalue], [int width], [char* format],
!      [log waves], [int lbound], [char* options] )

! Note that most of the optional parameters have default values
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
!       B              show Bandwidth, % of array that is non-zero
!       H              show rank, TheShape of array
!       L              laconic; skip printing name, size of array
!       N              Show where NaNs and Infs are located
!       R              rms       -- min, max, etc.
!       b              table of % vs. amount of differences (pdf)
!       c              clean
!       g              gaps      
!       l              collapse (last index)
!       r              ratios    -- min, max, etc. of differences' ratios
!       s              stats     -- number of differences
!       p              transpose 
!       t              trim      
!       u              unique    
!       w              wholearray
!       W[i]           wholearray, looping over ith index (for rank 3 and 4 arrays only)
!       1 or 2 or ..   ignored; calling routine is free to interpret

! An exception is the behavior of wholearray:
! if all {HRblrs} are FALSE, i.e. unset, the whole array is dumped (or diffed)
! if any is TRUE the whole array will be dumped only if
! w or wholearray is set to TRUE

! (for diff_fun)
!   character      meaning
!    ---           -------
!     a            diff the absolute values
!     r            divide the difference by the max abs of the 2 values
!     f            treat auxvalue as a fillvalue: return fillvalue if either
!                    value is fillvalue
!     p            treat auxvalue as a period: return 
!                       min(value1 - value2 + n*period)

! in the above, a string list is a string of elements (usu. comma-separated)
! === (end of api) ===

  public :: DIFF, DIFF_FUN, &
    & DUMP, DUMP_2x2xN, DUMPDATES, DUMPDUMPOPTIONS, DUMPLISTS, DUMPNAMEDVALUES, &
    & DUMPSUMS, DUMPTABLE, RESTOREDUMPCONFIG, SELFDIFF

  interface DIFF        ! dump diffs between pair of n-d arrays of numeric type
    module procedure DIFF_1D_DOUBLE, DIFF_1D_INTEGER, DIFF_1D_REAL
    module procedure DIFF_2D_DOUBLE, DIFF_2D_INTEGER, DIFF_2D_REAL
    module procedure DIFF_3D_DOUBLE, DIFF_3D_REAL
    module procedure DIFF_4D_DOUBLE, DIFF_4D_REAL
  end interface

  interface DIFF_FUN    ! return diffs between args or arrays of numeric type
    module procedure DIFF_SCALAR_DOUBLE, DIFF_SCALAR_REAL
  end interface

  interface FILTEREDDIFF        ! dump FILTEREDDIFFs between pair of n-d arrays of numeric type
    module procedure FILTEREDDIFF_1D_DOUBLE, FILTEREDDIFF_1D_INTEGER, FILTEREDDIFF_1D_REAL
    module procedure FILTEREDDIFF_2D_DOUBLE, FILTEREDDIFF_2D_INTEGER, FILTEREDDIFF_2D_REAL
    module procedure FILTEREDDIFF_3D_DOUBLE, FILTEREDDIFF_3D_REAL
    module procedure FILTEREDDIFF_4D_DOUBLE, FILTEREDDIFF_4D_REAL
  end interface

  interface DUMP        ! dump n-d arrays of homogeneous type
    module procedure DUMP_1D_BIT, DUMP_1D_CHAR, DUMP_1D_COMPLEX, DUMP_1D_DCOMPLEX
    module procedure DUMP_1D_DOUBLE, DUMP_1D_INTEGER, DUMP_1D_INTEGER_2B
    module procedure DUMP_1D_LOGICAL, DUMP_1D_REAL
    module procedure DUMP_2D_CHAR, DUMP_2D_COMPLEX, DUMP_2D_DCOMPLEX
    module procedure DUMP_2D_DOUBLE, DUMP_2D_INTEGER, DUMP_2D_INTEGER_2B
    module procedure DUMP_2D_LOGICAL, DUMP_2D_REAL
    module procedure DUMP_3D_CHAR, DUMP_3D_DOUBLE, DUMP_3D_INTEGER
    module procedure DUMP_3D_REAL, DUMP_3D_COMPLEX, DUMP_3D_DCOMPLEX
    module procedure DUMP_HASH_LOG, DUMP_HASH_STR, DUMP_STRLIST
    module procedure DUMP_4D_DOUBLE, DUMP_4D_REAL
  end interface

  interface DUMP_2x2xN ! For polarized incremental optical depth
    module procedure DUMP_2x2xN_COMPLEX, DUMP_2x2xN_DCOMPLEX
  end interface

  interface DUMPCOLLAPSEDARRAY
    module procedure DUMPCOLLAPSEDARRAY_1D_DOUBLE, DUMPCOLLAPSEDARRAY_1D_REAL
    module procedure DUMPCOLLAPSEDARRAY_2D_DOUBLE, DUMPCOLLAPSEDARRAY_2D_REAL
    module procedure DUMPCOLLAPSEDARRAY_3D_DOUBLE, DUMPCOLLAPSEDARRAY_3D_REAL
    module procedure DUMPCOLLAPSEDARRAY_1D_INTEGER
    module procedure DUMPCOLLAPSEDARRAY_2D_INTEGER
    module procedure DUMPCOLLAPSEDARRAY_3D_INTEGER
  end interface

  interface dumpDates
    module procedure dump_tai
  end interface

  interface DUMPLISTS
    module procedure DUMPLISTS_CHARS, DUMPLISTS_INTS
  end interface

  interface DUMPNAMEDVALUES   ! dump name-value pairs, names in string list
    module procedure DUMPNAMEDVALUES_DOUBLE, DUMPNAMEDVALUES_INTEGER
    module procedure DUMPNAMEDVALUES_REAL
  end interface

  interface DUMPTABLE   ! dump table of values, headers
    module procedure DUMPTABLE_DOUBLE, DUMPTABLE_INTEGER
    module procedure DUMPTABLE_REAL
  end interface

  interface SELFDIFF       ! dump increments between successive array values
    module procedure SELFDIFF_INTEGER
    module procedure SELFDIFF_REAL
    module procedure SELFDIFF_DOUBLE
  end interface

  interface DUMPSUMS       ! dump after summing successive array values
    module procedure DUMPSUMS_INTEGER
    module procedure DUMPSUMS_REAL
    module procedure DUMPSUMS_DOUBLE
  end interface

  interface PRINTIT
    module procedure PRINTIT_CHAR, PRINTIT_DOUBLE, PRINTIT_INT, PRINTIT_REAL
    module procedure PRINTIT_COMPLEX, PRINTIT_DCOMPLEX
  end interface

  interface PRINTRMSETC
    module procedure PRINTRMSETC_DOUBLE, PRINTRMSETC_INT, PRINTRMSETC_REAL
  end interface

  interface SAY_FILL
    module procedure SAY_FILL_CHAR, SAY_FILL_DOUBLE, SAY_FILL_INT
    module procedure SAY_FILL_REAL, SAY_FILL_COMPLEX, SAY_FILL_DCOMPLEX
  end interface

  interface UNFILTEREDDIFF        ! dump UNFILTEREDDIFFs between pair of n-d arrays of numeric type
    module procedure UNFILTEREDDIFF_1D_DOUBLE, UNFILTEREDDIFF_1D_INTEGER, UNFILTEREDDIFF_1D_REAL
    module procedure UNFILTEREDDIFF_2D_DOUBLE, UNFILTEREDDIFF_2D_INTEGER, UNFILTEREDDIFF_2D_REAL
    module procedure UNFILTEREDDIFF_3D_DOUBLE, UNFILTEREDDIFF_3D_REAL
    module procedure UNFILTEREDDIFF_4D_DOUBLE, UNFILTEREDDIFF_4D_REAL
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! These public parameters can't be reconfigured outside the module
  ! --------------------------------------------------------------------------
  character, public, parameter :: AFTERSUB = '#'
  character(len=*), parameter :: DEFAULTPCTFORMAT = '(0pf6.1)'
  ! These are the possible options to dumps, diffs
  character, public, parameter :: dopt_bandwidth   = 'B'
  character, public, parameter :: dopt_clean       = 'c'
  character, public, parameter :: dopt_collapse    = 'l'
  character, public, parameter :: dopt_cyclic      = 'y'
  character, public, parameter :: dopt_gaps        = 'g'
  character, public, parameter :: dopt_laconic     = 'L'
  character, public, parameter :: dopt_NaNs        = 'N'
  character, public, parameter :: dopt_ratios      = 'r'
  character, public, parameter :: dopt_rms         = 'R'
  character, public, parameter :: dopt_shape       = 'H'
  character, public, parameter :: dopt_stats       = 's'
  character, public, parameter :: dopt_table       = 'b'
  character, public, parameter :: dopt_transpose   = 'p'
  character, public, parameter :: dopt_trim        = 't'
  character, public, parameter :: dopt_unique      = 'u'
  character, public, parameter :: dopt_wholearray  = 'w'
  ! The following character strings can include one or more options listed above
  !
  ! E.g., '-crt' turns on CLEAN, RMS, and TRIM
  character(len=8), public, save :: DEFAULTDIFFOPTIONS = ' '
  character(len=8), public, save :: DEFAULTDUMPOPTIONS = ' '
  integer, public, save          :: DEFAULTMAXLON = 128
  character(len=8), public, save :: DUMPTABLESIDE      = 'top'
  logical, public, save ::   DIFFRMSMEANSRMS           = .false.
  logical, public, save ::   DONTDUMPIFALLEQUAL        = .true.
  logical, public, save ::   FILTERFILLSFROMRMS        = .false.
  logical, public, save ::   PRINTFILLVALUE            = .true.
  logical, public, save ::   PRINTNAMEIFDIFF           = .true.
  logical, public, save ::   STATSONONELINE            = .true.

  ! This determines how a higher-rank array is collapsed to a lower-rank one
  character(len=16), public, save :: COLLAPSEOPTIONS = 'num[+]all[+]'

  ! These determine how dumped numerical data (s.p. or d.p.) will be formatted
  character(len=2), public, save  :: INTPLACES = '6' ! how many places
  integer,          public, save  :: MAXNUMNANS= 60  ! how many NaNs to show
  character(len=16), public, save :: PCTFORMAT = '*' ! * means default format
  character(len=16), public, save :: RMSFORMAT = '*' ! * means default format
  character(len=16), public, save :: SDFORMATDEFAULT = '(1pg14.6)'
  character(*), parameter :: sdFormatDefaultCmplx = &
    & '(1x,"(",1pg13.6,",",1pg13.6,")")'
  ! --------------------------------------------------------------------------

  ! These are private variables declared module-wide purely for convenience
  integer, parameter :: MAXLINELEN = 120
  integer, parameter :: MAXNUMELEMENTS = 2000
  integer, parameter :: TOOMANYELEMENTS = 125*50*3500 ! Don't try to diff l1b DACS
  logical, parameter ::   DEEBUG = .false.
  logical, parameter ::   SHORTCUTDIFFS = .false.
  ! character(len=MAXLINELEN) :: LINEOFZEROS
  logical :: DUMPTHESEZEROS
  logical :: myBandwidth, myClean, myCollapse, myCyclic, myDirect, myGaps, &
    & myLaconic, myNaNs, myRatios, myRMS, myShape, myStats, &
    & myTable, myTranspose, myTrim, myUnique, myWholeArray, onlyWholeArray
  character(len=16) :: myOptions
  character(len=16) :: myPCTFormat
  logical, save :: nameHasBeenPrinted = .false.
  integer :: bwidth, myRank, numNonFill, numFill, indx2BSliced, iSlice
  real :: pctnzero
  logical, save :: thisIsADiff = .false.
  integer :: how_many
  integer, dimension(1024) :: which
  complex, parameter :: one_c4 = (1., 0.)

contains

 ! ---------------------------------------------  DIFF  -----
 ! This family of routines dumps the differences between two arrays
 ! with the same shape and numeric type
 ! Its behavior is modified by the following optional parameters
 ! FillValue   Ignore these values when computing min, max
 ! Width       Row size when dumping wholearrays
 ! Format      Output format when printing wholearray
 ! LBound      Lower bound when printing wholearray indices
 !
 ! The following optional params are now set by options
 ! Wholearray  Whether to print whole array of differences
 ! Stats       Dump number of differences found, %
 ! RMS         Dump min, max, rms
 ! Clean       Clean up after any prior dumps
  subroutine DIFF_1D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    double precision, intent(in) :: ARRAY1(:)
    character(len=*), intent(in) :: NAME1
    double precision, intent(in) :: ARRAY2(:)
    character(len=*), intent(in) :: NAME2
    double precision, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(FillValue) ) then
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    elseif ( product(shape(array1)) > TOOMANYELEMENTS ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    else
      call FilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    endif
  end subroutine DIFF_1D_DOUBLE

  subroutine DIFF_1D_INTEGER ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, intent(in) :: ARRAY1(:)
    character(len=*), intent(in) :: NAME1
    integer, intent(in) :: ARRAY2(:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(FillValue) ) then
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    elseif ( product(shape(array1)) > TOOMANYELEMENTS ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    else
      call FilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    endif
  end subroutine DIFF_1D_INTEGER

  subroutine DIFF_1D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    real, intent(in) :: ARRAY1(:)
    character(len=*), intent(in) :: NAME1
    real, intent(in) :: ARRAY2(:)
    character(len=*), intent(in) :: NAME2
    real, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(FillValue) ) then
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    elseif ( product(shape(array1)) > TOOMANYELEMENTS ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    else
      call FilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    endif
  end subroutine DIFF_1D_REAL

  subroutine DIFF_2D_INTEGER ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, intent(in) :: ARRAY1(:,:)
    character(len=*), intent(in) :: NAME1
    integer, intent(in) :: ARRAY2(:,:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(FillValue) ) then
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    elseif ( product(shape(array1)) > TOOMANYELEMENTS ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    else
      call FilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    endif
  end subroutine DIFF_2D_INTEGER

  subroutine DIFF_2D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    double precision, intent(in) :: ARRAY1(:,:)
    character(len=*), intent(in) :: NAME1
    double precision, intent(in) :: ARRAY2(:,:)
    character(len=*), intent(in) :: NAME2
    double precision, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(FillValue) ) then
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    elseif ( product(shape(array1)) > TOOMANYELEMENTS ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    else
      call FilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    endif
  end subroutine DIFF_2D_DOUBLE

  subroutine DIFF_2D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    real, intent(in) :: ARRAY1(:,:)
    character(len=*), intent(in) :: NAME1
    real, intent(in) :: ARRAY2(:,:)
    character(len=*), intent(in) :: NAME2
    real, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(FillValue) ) then
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    elseif ( product(shape(array1)) > TOOMANYELEMENTS ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    else
      call FilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    endif
  end subroutine DIFF_2D_REAL

  subroutine DIFF_3D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    double precision, intent(in) :: ARRAY1(:,:,:)
    character(len=*), intent(in) :: NAME1
    double precision, intent(in) :: ARRAY2(:,:,:)
    character(len=*), intent(in) :: NAME2
    double precision, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(FillValue) ) then
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    elseif ( product(shape(array1)) > TOOMANYELEMENTS ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    else
      call FilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    endif
  end subroutine DIFF_3D_DOUBLE

  subroutine DIFF_3D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    real, intent(in) :: ARRAY1(:,:,:)
    character(len=*), intent(in) :: NAME1
    real, intent(in) :: ARRAY2(:,:,:)
    character(len=*), intent(in) :: NAME2
    real, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(FillValue) ) then
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    elseif ( product(shape(array1)) > TOOMANYELEMENTS ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    else
      call FilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    endif
  end subroutine DIFF_3D_REAL

  subroutine DIFF_4D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    double precision, intent(in) :: ARRAY1(:,:,:,:)
    character(len=*), intent(in) :: NAME1
    double precision, intent(in) :: ARRAY2(:,:,:,:)
    character(len=*), intent(in) :: NAME2
    double precision, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    call theDumpBegins ( options )
    if ( any(shape(array1) == 0) ) then
      call output( 'array sizes are 0', advance='yes' )
    else if ( size(array1,1) == 1 ) then
      if ( myWholeArray ) call output( '1st size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(1,:,:,:), name1, array2(1,:,:,:), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( size(array1,2) == 1 ) then
      if ( myWholeArray ) call output( '2nd size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(:,1,:,:), name1, array2(:,1,:,:), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( size(array1,3) == 1 ) then
      if ( myWholeArray ) call output( '3rd size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(:,:,1,:), name1, array2(:,:,1,:), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( size(array1,4) == 1 ) then
      if ( myWholeArray ) call output( '4th size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(:,:,:,1), name1, array2(:,:,:,1), name2, &
        & fillvalue, width, format, lbound, options )
    elseif ( .not. present(FillValue) ) then
      if ( myWholeArray ) call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    elseif ( product(shape(array1)) > TOOMANYELEMENTS ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    else
      call FilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    endif
  end subroutine DIFF_4D_DOUBLE

  subroutine DIFF_4D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    real, intent(in) :: ARRAY1(:,:,:,:)
    character(len=*), intent(in) :: NAME1
    real, intent(in) :: ARRAY2(:,:,:,:)
    character(len=*), intent(in) :: NAME2
    real, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    call theDumpBegins ( options )
    if ( any(shape(array1) == 0) ) then
      call output( 'array sizes are 0', advance='yes' )
    else if ( size(array1,1) == 1 ) then
      if ( myWholeArray ) call output( '1st size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(1,:,:,:), name1, array2(1,:,:,:), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( size(array1,2) == 1 ) then
      if ( myWholeArray ) call output( '2nd size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(:,1,:,:), name1, array2(:,1,:,:), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( size(array1,3) == 1 ) then
      if ( myWholeArray ) call output( '3rd size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(:,:,1,:), name1, array2(:,:,1,:), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( size(array1,4) == 1 ) then
      if ( myWholeArray ) call output( '4th size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(:,:,:,1), name1, array2(:,:,:,1), name2, &
        & fillvalue, width, format, lbound, options )
    elseif ( .not. present(FillValue) ) then
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    elseif ( product(shape(array1)) > TOOMANYELEMENTS ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & WIDTH, FORMAT, LBOUND, OPTIONS )
    else
      call FilteredDiff( ARRAY1, NAME1, ARRAY2, NAME2, &
        & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    endif
  end subroutine DIFF_4D_REAL

  ! -----------------------------------------------  DIFF_SCALAR  -----
  ! This family of functions differences two values and returns
  ! according to options:
  ! if options contains
  ! character        meaning
  ! ---------        -------
  !     a            diff the absolute values
  !     r            divide the difference by the max abs of the 2 values
  !     f            treat auxvalue as a fillvalue: return fillvalue if either
  !                    value is fillvalue
  !     p            treat auxvalue as a period: return 
  !                    min(value1 - value2 + n*period)
  elemental function DIFF_SCALAR_DOUBLE ( VALUE1, VALUE2, AUXVALUE, OPTIONS ) &
    & result(d)
    ! Args
    double precision, intent(in)           :: value1
    double precision, intent(in)           :: value2
    double precision, optional, intent(in) :: auxvalue
    character(len=*), optional, intent(in) :: options
    double precision                       :: d
    ! Internal variables
    integer, parameter :: RK = kind(0.0d0)
    ! Executable
    include 'diff_scalar.f9h'
  end function DIFF_SCALAR_DOUBLE

  elemental function DIFF_SCALAR_REAL ( VALUE1, VALUE2, AUXVALUE, OPTIONS ) &
    & result(d)
    ! Args
    real, intent(in)                       :: value1
    real, intent(in)                       :: value2
    real, optional, intent(in)             :: auxvalue
    character(len=*), optional, intent(in) :: options
    real                                   :: d
    ! Internal variables
    integer, parameter :: RK = kind(0.0e0)
    ! Executable
    include 'diff_scalar.f9h'
  end function DIFF_SCALAR_REAL

  ! -----------------------------------------------  DUMP_1D_BIT  -----
  subroutine DUMP_1D_BIT ( ARRAY, NAME, BITNAMES, FILLVALUE, OPTIONS )
    integer, intent(in) :: ARRAY(:)
    character(len=*), intent(in) :: NAME
    character(len=*), intent(in) :: BITNAMES
    integer, intent(in), optional :: FILLVALUE
    character(len=*), optional, intent(in) :: options

    integer :: howMany, J, K
    integer, dimension(MAXBITNUMBER) :: ints
    integer :: myFillValue
    integer :: NumBitNames
    integer :: NumZeroRows
    integer, dimension(MAXBITNUMBER+1) :: set
    ! Executable
    call theDumpBegins ( options )
    myFillValue = 0
    if ( present(FillValue) ) myFillValue = FillValue

    NumBitNames = NumStringElements( bitNames, countEmpty=.true. )
    if ( NumBitNames < 1 ) NumBitNames = MAXBITNUMBER+1
    NumBitNames = min( numBitNames, MAXBITNUMBER+1 )
    numZeroRows = 0
    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else
      call name_and_size ( name, myClean, size(array) )
      call newline
      call output( trim(BitNames), advance='yes' )
      do j=1, size(array)
        if ( array(j) == myFillValue ) then
          numZeroRows = numZeroRows + 1
        else
          call WhichBitsAreSet( array(j), set, howMany )
          ints = 0
          do k=1, howMany
            ints( 1 + set(k) ) = 1
          enddo
          if ( numZeroRows > 0 ) then
            call output ( ' ' , advance='no' )
            call output ( numZeroRows , advance='no' )
            call output ( ' lines of ', advance='no' )
            call output ( myFillValue , advance='no' )
            call output ( ' not printed', advance='yes' )
          endif
          call output( ints(1:numBitNames), format='(i3)', advance='yes' )
          numZeroRows = 0
        endif
      enddo
    endif
    call theDumpEnds
  end subroutine DUMP_1D_BIT

  ! -----------------------------------------------  DUMP_1D_CHAR  -----
  subroutine DUMP_1D_CHAR ( ARRAY, NAME, FILLVALUE, WIDTH, OPTIONS, MAXLON, &
    & TheSHAPE )
    character(len=*), intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), optional, intent(in) :: OPTIONS
    integer, intent(in), optional :: MAXLON
    character(len=*), intent(in), optional :: TheShape

    integer :: J, K
    integer :: LON
    integer :: NumZeroRows
    character(len=len(array)) :: myFillValue
    integer :: MyWidth

    call theDumpBegins ( options )
    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue
    MyWidth = 10
    if ( present(width) ) MyWidth = width

    lon = len(array(1))
    if ( myTrim ) lon = maxval(len_trim(array))
    if ( present(maxlon) ) then
      lon = min(lon, maxlon)
    else if ( .not. myTrim ) then
      lon = min(lon, DEFAULTMAXLON)
    end if

    numZeroRows = 0
    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1, TheShape )
      call output ( array(1)(1:lon), advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) .and. .not. mylaconic ) call newLine
      do j = 1, size(array), MyWidth
        DumpTheseZeros = myClean .or. &
          & any(array(j:min(j+2*myWidth-1, size(array))) /= myFillValue)
        if (.not. myClean) then
          if ( DumpTheseZeros ) then
            call say_fill ( (/ j-1, size(array) /), numZeroRows, &
              & myFillValue, inc=1 )
          else
            numZeroRows = numZeroRows + 1
          end if
        end if
        if ( DumpTheseZeros ) then
          do k = j, min(j+MyWidth-1, size(array))
              call output ( array(k)(1:lon) // ' ' , advance='no' )
          end do
          call newLine
          numZeroRows = 0
        end if
      end do ! j
      call say_fill ( (/ j-MyWidth, size(array) /), numZeroRows, &
        & myFillValue )
    end if
    call theDumpEnds
  end subroutine DUMP_1D_CHAR

  ! --------------------------------------------  DUMP_1D_COMPLEX  -----
  subroutine DUMP_1D_COMPLEX ( ARRAY, NAME, WIDTH, FORMAT, &
    & FILLVALUE, LBOUND, OPTIONS, TheShape )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    real(rk), intent(in), optional :: FillValue
    integer, intent(in), optional :: LBOUND ! Low bound for Array             
    character(len=*), optional, intent(in) :: OPTIONS
    character(len=*), intent(in), optional :: TheShape

    integer :: J, K, MyWidth
    character(len=64) :: MyFormat
    complex(rk) :: myFillValue
    integer :: BASE
    integer :: NumZeroRows
    ! Executable
    myFormat = sdFormatDefaultCmplx
    call theDumpBegins ( options )
    myFillValue = 0.
    myWidth = 3
    if ( present(width) ) myWidth = width
    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )
    base = 0
    if ( present(lbound) ) base = lbound - 1

    include 'dump1db.f9h'
  end subroutine DUMP_1D_COMPLEX

  ! --------------------------------------------  DUMP_1D_DCOMPLEX  -----
  subroutine DUMP_1D_DCOMPLEX ( ARRAY, NAME, WIDTH, FORMAT, &
    & FILLVALUE, LBOUND, OPTIONS, TheShape )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    real(rk), intent(in), optional :: FillValue
    integer, intent(in), optional :: LBOUND ! Low bound for Array             
    character(len=*), optional, intent(in) :: OPTIONS
    character(len=*), intent(in), optional :: TheShape

    integer :: J, K, MyWidth
    character(len=64) :: MyFormat
    complex(rk) :: myFillValue
    integer :: BASE
    integer :: NumZeroRows
    ! Executable
    myFormat = sdFormatDefaultCmplx
    call theDumpBegins ( options )
    myFillValue = 0.
    myWidth = 3
    if ( present(width) ) myWidth = width
    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )
    base = 0
    if ( present(lbound) ) base = lbound - 1

    include 'dump1db.f9h'
  end subroutine DUMP_1D_DCOMPLEX

 ! ---------------------------------------------  DUMP_1D_DOUBLE  -----
  subroutine DUMP_1D_DOUBLE ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS, TheShape )
    double precision, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), optional, intent(in) :: OPTIONS
    character(len=*), intent(in), optional :: TheShape

    integer, dimension(MAXNUMELEMENTS) :: counts
    double precision, dimension(MAXNUMELEMENTS) :: elements
    double precision :: myFillValue
    integer :: Base, J, K
    character(len=64) :: MyFormat
    integer :: MyWidth
    integer :: NumZeroRows
    integer :: nUnique
    myFormat = sdFormatDefault
    include 'dump1d.f9h'
    include 'dump1db.f9h'
  end subroutine DUMP_1D_DOUBLE

  ! --------------------------------------------  DUMP_BOGUS  -----
  ! Never used--just here to tell Makefiles that dumpstats.f9h is
  ! a prerequisite for dump_0 because perl script f90makedep.pl
  ! won't follow .f9h files to look for more uses and includes.
  ! When will we repair the perl script?
  !
  subroutine DUMP_BOGUS ( ARRAY, NAME, &
    & FILLVALUE, FORMAT, WIDTH, LBOUND, OPTIONS )
    integer, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many numbers per line (10)?
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), optional, intent(in) :: options
    integer :: myFillValue
    myFillValue = 0
    myRank = 0
    include 'dumpstats.f9h'
  end subroutine DUMP_BOGUS

  ! --------------------------------------------  DUMP_1D_INTEGER  -----
  subroutine DUMP_1D_INTEGER ( ARRAY, NAME, &
    & FILLVALUE, FORMAT, WIDTH, LBOUND, OPTIONS, TheShape )
    integer, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many numbers per line (10)?
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), optional, intent(in) :: OPTIONS
    character(len=*), intent(in), optional :: TheShape

    integer, dimension(MAXNUMELEMENTS) :: counts
    integer, dimension(MAXNUMELEMENTS) :: elements
    integer :: myFillValue
    integer :: Base, J, K
    character(len=64) :: MyFormat
    integer :: MyWidth
    integer :: NumZeroRows
    integer :: nUnique
    myFormat = 'places=' // INTPLACES ! To sneak places arg into call to output
    include 'dump1d.f9h'
    include 'dump1db.f9h'
  end subroutine DUMP_1D_INTEGER

  ! --------------------------------------------  DUMP_1D_INTEGER_2B  -----
  subroutine DUMP_1D_INTEGER_2B ( ARRAY, NAME, &
    & FILLVALUE, FORMAT, WIDTH, LBOUND, OPTIONS, TheShape )
    use ISO_C_BINDING, only: C_int16_t
    integer(C_int16_t), intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many numbers per line (10)?
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), optional, intent(in) :: OPTIONS
    character(len=*), intent(in), optional :: TheShape

    integer, dimension(MAXNUMELEMENTS) :: counts
    integer, dimension(MAXNUMELEMENTS) :: elements
    integer :: myFillValue
    integer :: Base, J, K
    character(len=64) :: MyFormat
    integer :: MyWidth
    integer :: NumZeroRows
    integer :: nUnique
    call dump( int(array), NAME, &
    & FILLVALUE, FORMAT, WIDTH, LBOUND, OPTIONS, TheShape )
  end subroutine DUMP_1D_INTEGER_2B

  ! ----------------------------------------------  DUMP_1D_LOGICAL ----
  subroutine DUMP_1D_LOGICAL ( ARRAY, NAME, LBOUND, OPTIONS, TheShape )
    logical, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: LBOUND ! Low bound of Array
    character(len=*), optional, intent(in) :: OPTIONS
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, J, K, N, NTRUE, NFALSE
    ! Executable

    call theDumpBegins ( options )
    base = 0
    if ( present(lbound) ) base = lbound - 1

    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else if ( size(array) == 1 .and. base == 0 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), advance='yes' )
    else if ( all(array .eqv. .true.) ) then
      call name_and_size ( name, myClean, size(array) )
      call output ( ': all ', advance='no' )
      call output( trim(arrayShapeToString(shape(array))), advance='no' )
      call output( ' values are T', advance='yes' )
    else if ( all(array .eqv. .false.) ) then
      call name_and_size ( name, myClean, size(array) )
      call output ( ': all ', advance='no' )
      call output( trim(arrayShapeToString(shape(array))), advance='no' )
      call output( ' values are F', advance='yes' )
    elseif ( myGaps ) then
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) .and. .not. mylaconic ) call newLine
      k = 0
      if ( size(array)/100 > 100 )  call showColumnNums( 100 )
      do j=1, size(array), 100
        N = min(k+100, size(array)) - k
        call output( k+1 , advance='no' )
        call output ( ' through ', advance='no' )
        call output( k+N, advance='yes' )
        ntrue = count( array(k+1 : min(k+100, size(array))) )
        nfalse = n - ntrue
        if ( ntrue > 0 .and. ( ntrue < nfalse .or. nfalse < 1 ) ) then
          call output( array(k+1 : k+N), onlyif=.true., advance='yes' )
        else
          call output( array(k+1 : k+N), onlyif=.false., advance='yes' )
        endif
        k = k + 100
      enddo
      call showColumnNums( 100 )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) .and. .not. mylaconic ) call newLine
      do j = 1, size(array), 34
        if (.not. myClean) then
          call output ( j+base, max(4,ilog10(size(array))+1) , advance='no' )
          call output ( afterSub , advance='no' )
        end if
        do k = j, min(j+33, size(array))
          call output ( array(k) , advance='no' )
        end do
        call newLine
      end do
    end if
    call theDumpEnds
  end subroutine DUMP_1D_LOGICAL

  ! -----------------------------------------------  DUMP_1D_REAL  -----
  subroutine DUMP_1D_REAL ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS, TheShape )
    real, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    real, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), optional, intent(in) :: OPTIONS
    character(len=*), intent(in), optional :: TheShape

    integer, dimension(MAXNUMELEMENTS) :: counts
    real, dimension(MAXNUMELEMENTS) :: elements
    real :: myFillValue
    integer :: Base, J, K
    character(len=64) :: MyFormat
    integer :: MyWidth
    integer :: NumZeroRows
    integer :: nUnique
    myFormat = sdFormatDefault
    include 'dump1d.f9h'
    include 'dump1db.f9h'
  end subroutine DUMP_1D_REAL

  ! -----------------------------------------------  DUMP_2D_CHAR  -----
  recursive subroutine DUMP_2D_CHAR ( ARRAY, NAME, FILLVALUE, WIDTH, MAXLON, &
    & OPTIONS, TheShape )
    character(len=*), intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: MAXLON
    character(len=*), optional, intent(in) :: options
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: TheShape

    integer :: I, J, K
    integer :: LON
    integer :: NumZeroRows
    character(len=len(array)) :: myFillValue
    integer :: MyWidth

    call theDumpBegins ( options )
    if ( myTranspose ) then
      call dump ( transpose(array), name, &
        & fillvalue, width, &
        & options=snipoption(options, dopt_transpose) )
      return
    endif
    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue
    MyWidth = 10
    if ( present(width) ) MyWidth = width
    lon = len(array(1,1))
    if ( myTrim ) lon = maxval(len_trim(array))
    if ( present(maxlon) ) lon = min(lon, maxlon)

    numZeroRows = 0
    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else if ( size(array,1) == 1 ) then
      call name_and_size ( name, myClean, 1, TheShape )
      call output ( array(1,:), advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, fillValue=fillValue, maxlon=maxlon, options=options )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) .and. .not. mylaconic ) call newLine
      do i = 1, size(array,1)
        do j = 1, size(array,2), MyWidth
          DumpTheseZeros = myClean .or. &
            & any(array(i,j:min(j+2*MyWidth-1, size(array,2))) /= myFillValue) &
            & .or. &
            & ( j+MyWidth >= size(array,2) .and. &
            & any(array(min(i+1, size(array,1)),1:min(1+MyWidth-1, size(array,2))) &
            & /= myFillValue) &
            & )
          if (.not. myClean) then
            if ( DumpTheseZeros ) then
              call say_fill ( (/ i, size(array,1), j-1, size(array,2) /), &
                & numZeroRows, myFillValue, inc=3 )
            else
              numZeroRows = numZeroRows + 1
            end if
          end if
          if ( DumpTheseZeros ) then
            do k = j, min(j+myWidth-1, size(array,2))
                call output ( array(i,k)(1:lon) // ' ' , advance='no' )
            end do
            call newLine
            numZeroRows = 0
          end if
        end do ! j
      end do ! i
      call say_fill ( (/ i-1, size(array,1), j-MyWidth, size(array,2) /), &
        & numZeroRows, myFillValue )
    end if
    call theDumpEnds
  end subroutine DUMP_2D_CHAR

  ! --------------------------------------------  DUMP_2D_COMPLEX  -----
  recursive subroutine DUMP_2D_COMPLEX ( ARRAY, NAME, WIDTH, FORMAT, &
    & FILLVALUE, OPTIONS, LBOUND, TheShape )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: WIDTH ! How many per line?
    character(len=*), optional :: FORMAT
    real, intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: options
    integer, intent(in), optional :: LBound ! to print for first dimension
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, J, K
    integer :: myWidth
    integer :: NumZeroRows
    complex(rk) :: myFillValue
    character(len=64) :: MyFormat
    ! Executable

    myFormat = sdFormatDefaultCmplx
    call theDumpBegins ( options )
    if ( myTranspose ) then
      call dump ( transpose(array), name, &
        & width, format, fillvalue, &
        & options=snipoption(options, dopt_transpose) )
      return
    endif

    myWidth = 3
    if ( present(width) ) myWidth = width

    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )

    myFillValue = 0.0_rk
    if ( present(fillValue) ) myFillValue = fillValue

    base = 1
    if ( present(lbound) ) base = lbound

    include 'dump2db.f9h'
  end subroutine DUMP_2D_COMPLEX

  ! --------------------------------------------  DUMP_2D_DCOMPLEX  -----
  recursive subroutine DUMP_2D_DCOMPLEX ( ARRAY, NAME, WIDTH, FORMAT, &
    & FILLVALUE, OPTIONS, LBOUND, TheShape )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: WIDTH ! How many per line?
    character(len=*), optional :: FORMAT
    real(rk), intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: options
    integer, intent(in), optional :: LBound ! to print for first dimension
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, J, K
    integer :: myWidth
    integer :: NumZeroRows
    complex(rk) :: myFillValue
    character(len=64) :: MyFormat

    ! Executable

    myFormat = sdFormatDefaultCmplx
    call theDumpBegins ( options )
    if ( myTranspose ) then
      call dump ( transpose(array), name, &
        & width, format, fillvalue, &
        & options=snipoption(options, dopt_transpose) )
      return
    endif

    myWidth = 3
    if ( present(width) ) myWidth = width

    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )

    myFillValue = 0.0_rk
    if ( present(fillValue) ) myFillValue = fillValue

    base = 1
    if ( present(lbound) ) base = lbound

    include 'dump2db.f9h'
  end subroutine DUMP_2D_DCOMPLEX

  ! ---------------------------------------------  DUMP_2D_DOUBLE  -----
 recursive subroutine DUMP_2D_DOUBLE ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS, TheShape )
    double precision, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! to print for first dimension
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, J, K
    integer :: NumZeroRows
    double precision :: myFillValue
    character(len=64) :: MyFormat
    integer :: nUnique
    integer :: MyWidth
    integer, dimension(MAXNUMELEMENTS) :: counts
    double precision, dimension(MAXNUMELEMENTS) :: elements

    myFormat = sdFormatDefault

    base = 1
    if ( present(lbound) ) base = lbound

    include 'dump2d.f9h'
    include 'dump2db.f9h'
  end subroutine DUMP_2D_DOUBLE

  ! --------------------------------------------  DUMP_2D_INTEGER  -----
  recursive subroutine DUMP_2D_INTEGER ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS, TheShape )
    integer, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many numbers per line (10)?
    integer, intent(in), optional :: LBOUND ! to print for first dimension
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, J, K
    integer :: MyWidth
    integer :: NumZeroRows
    integer :: myFillValue
    character(len=64) :: MyFormat
    integer :: nUnique
    integer, dimension(MAXNUMELEMENTS) :: counts
    integer, dimension(MAXNUMELEMENTS) :: elements

    myFormat = 'places=' // INTPLACES ! To sneak places arg into call to output

    base = 1
    if ( present(lbound) ) base = lbound

    include 'dump2d.f9h'
    include 'dump2db.f9h'
  end subroutine DUMP_2D_INTEGER

  ! --------------------------------------------  DUMP_2D_INTEGER_2B  -----
  recursive subroutine DUMP_2D_INTEGER_2B ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS, TheShape )
    use ISO_C_BINDING, only: C_int16_t
    integer(C_int16_t), intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many numbers per line (10)?
    integer, intent(in), optional :: LBOUND ! to print for first dimension
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, J, K
    integer :: MyWidth
    integer :: NumZeroRows
    integer :: myFillValue
    character(len=64) :: MyFormat
    integer :: nUnique
    integer, dimension(MAXNUMELEMENTS) :: counts
    integer, dimension(MAXNUMELEMENTS) :: elements

    call dump( int(array), NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS, TheShape )
  end subroutine DUMP_2D_INTEGER_2B

  ! --------------------------------------------  DUMP_2D_LOGICAL  -----
  recursive subroutine DUMP_2D_LOGICAL ( ARRAY, NAME, OPTIONS, THESHAPE )
    logical, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: I, J, K
    integer, parameter :: MyWidth = 34

    call theDumpBegins ( options )
    if ( myTranspose ) then
      call dump ( transpose(array), name, &
        & options=snipoption(options, dopt_transpose) )
      return
    endif

    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1, TheShape )
      call output ( array(1,1), advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, options=options )
    else if ( all(array .eqv. .true.) ) then
      call name_and_size ( name, myClean, size(array), TheShape )
      call output ( ': all ', advance='no' )
      call output( trim(arrayShapeToString(shape(array))), advance='no' )
      call output( ' values are T', advance='yes' )
    else if ( all(array .eqv. .false.) ) then
      call name_and_size ( name, myClean, size(array), TheShape )
      call output ( ': all ', advance='no' )
      call output( trim(arrayShapeToString(shape(array))), advance='no' )
      call output( ' values are F', advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) .and. .not. mylaconic ) call newLine
      do i = 1, size(array,1)
        do j = 1, size(array,2), myWidth
          if (.not. myClean) then
            call output ( i, places=max(4,ilog10(size(array,1))+1) , advance='no' )
            call output ( j, places=max(4,ilog10(size(array,2))+1) , advance='no' )
            call output ( afterSub , advance='no' )
          end if
          do k = j, min(j+myWidth-1, size(array,2))
            call output ( array(i,k) , advance='no' )
          end do
          call newLine
        end do ! j
      end do ! i
    end if
    call theDumpEnds
  end subroutine DUMP_2D_LOGICAL

  ! -----------------------------------------------  DUMP_2D_REAL  -----
  recursive subroutine DUMP_2D_REAL ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS, THESHAPE )
    real, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    real, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, J, K
    integer :: NumZeroRows
    real :: myFillValue
    character(len=64) :: MyFormat
    integer :: nUnique
    integer :: MyWidth
    integer, dimension(MAXNUMELEMENTS) :: counts
    real, dimension(MAXNUMELEMENTS) :: elements

    myFormat = sdFormatDefault

    base = 1
    if ( present(lbound) ) base = lbound

    include 'dump2d.f9h'
    include 'dump2db.f9h'
  end subroutine DUMP_2D_REAL

  ! -----------------------------------------  DUMP_2x2xN_COMPLEX  -----
  subroutine DUMP_2x2xN_COMPLEX ( ARRAY, NAME, FORMAT, OPTIONS )
  ! This is for dumping polarized incremental optical depth
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: ARRAY(:,:,:) ! Better be 2x2xn
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: options

    integer :: J
    character(len=64) :: MyFormat

    call theDumpBegins ( options )
    myFormat = sdFormatDefaultCmplx
    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )

    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) .and. .not. mylaconic ) call newLine
      do j = 1, size(array,3)
        if (.not. myClean) then
          call output ( j, max(3,ilog10(size(array))+1) , advance='no' )
          call output ( afterSub , advance='no' )
        end if
        call output ( array(1,1,j), myFormat , advance='no' )
        call output ( array(1,2,j), myFormat, advance='yes' )
        call blanks ( max(3,ilog10(size(array))+1) + len(afterSub) )
        call output ( array(2,1,j), myFormat , advance='no' )
        call output ( array(2,2,j), myFormat, advance='yes' )
      end do
    end if
    call theDumpEnds
  end subroutine DUMP_2x2xN_COMPLEX

  ! ----------------------------------------  DUMP_2x2xN_DCOMPLEX  -----
  subroutine DUMP_2x2xN_DCOMPLEX ( ARRAY, NAME, FORMAT, OPTIONS )
  ! This is for dumping polarized incremental optical depth
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: ARRAY(:,:,:) ! Better be 2x2xn
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: options

    integer :: J
    character(len=64) :: MyFormat

    call theDumpBegins ( options )
    myFormat = sdFormatDefaultCmplx
    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )

    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) .and. .not. mylaconic ) call newLine
      do j = 1, size(array,3)
        if (.not. myClean) then
          call output ( j, max(3,ilog10(size(array))+1) , advance='no' )
          call output ( afterSub , advance='no' )
        end if
        call output ( array(1,1,j), myFormat , advance='no' )
        call output ( array(1,2,j), myFormat, advance='yes' )
        call blanks ( max(3,ilog10(size(array))+1) + len(afterSub) )
        call output ( array(2,1,j), myFormat , advance='no' )
        call output ( array(2,2,j), myFormat, advance='yes' )
      end do
    end if
    call theDumpEnds
  end subroutine DUMP_2x2xN_DCOMPLEX

  ! -----------------------------------------------  DUMP_3D_CHAR  -----
  subroutine DUMP_3D_CHAR ( ARRAY, NAME, FILLVALUE, WIDTH, &
    & MAXLON, OPTIONS, THESHAPE )
    character(len=*), intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    integer, intent(in), optional :: MAXLON
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: LON
    integer :: I, J, K, L
    integer :: NumZeroRows
    character(len=len(array)) :: myFillValue
    integer :: MyWidth

    call theDumpBegins ( options )
    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue

    lon = len(array(1,1,1))
    if ( myTrim ) lon = maxval(len_trim(array))
    if ( present(maxlon) ) lon = min(lon, maxlon)
    MyWidth = 10
    if ( present(width) ) MyWidth = width

    numZeroRows = 0
    if ( any(shape(array) == 0) ) then
      call empty ( name )
    else if ( size(array,1) == 1 ) then
      call dump ( array(1,:,:), name, fillValue=fillValue, &
        & maxlon=maxlon, options=options, width=width )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1,:), name, fillValue=fillValue, &
        & maxlon=maxlon, options=options, width=width )
    else if ( size(array,3) == 1 ) then
      call dump ( array(:,:,1), name, fillValue=fillValue, &
        & maxlon=maxlon, options=options, width=width )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) .and. .not. mylaconic ) call newLine
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), MyWidth
            DumpTheseZeros = myClean .or. &
              & any(array(i,j,k:min(k+2*MyWidth-1, size(array,3))) /= myFillValue) &
              & .or. &
              & ( k+MyWidth >= size(array,3) .and. &
              & any(array(i,min(j+1, size(array,2)),1:min(1+MyWidth-1, size(array,3))) &
              & /= myFillValue) &
              & )
            if (.not. myClean) then
              if ( DumpTheseZeros ) then
                call say_fill ( (/ i, size(array,1), j-1, size(array,2), &
                  & k, size(array,3) /), numZeroRows, myFillValue, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
          if ( DumpTheseZeros ) then
              do l = k, min(k+MyWidth-1, size(array,3))
                  call output ( array(i,j,l)(1:lon) // ' ' , advance='no' )
              end do
              call newLine
              numZeroRows = 0
            end if
          end do
        end do
      end do
      call say_fill ( (/ i-1, size(array,1), j-1, size(array,2), &
        & k-MyWidth, size(array,3) /), numZeroRows, myFillValue )
    end if
    call theDumpEnds
  end subroutine DUMP_3D_CHAR

  ! --------------------------------------------  DUMP_3D_COMPLEX  -----
  subroutine DUMP_3D_COMPLEX ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS, THESHAPE )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    real, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: THESHAPE

    integer :: Base, I, J, K, L
    integer :: NumZeroRows
    complex(rk) :: myFillValue
    character(len=64) :: MyFormat
    integer :: MyWidth

    ! Executable
    myFormat = sdFormatDefaultCmplx
    call theDumpBegins ( options )
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue

    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )
    
    MyWidth = 3
    if ( present(width) ) MyWidth = width

    include 'dump3db.f9h'
  end subroutine DUMP_3D_COMPLEX

  ! -------------------------------------------  DUMP_3D_DCOMPLEX  -----
  subroutine DUMP_3D_DCOMPLEX ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS, THESHAPE )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    real(rk), intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, J, K, L
    integer :: NumZeroRows
    complex(rk) :: myFillValue
    character(len=64) :: MyFormat
    integer :: myWidth

    ! Executable
    myFormat = sdFormatDefaultCmplx
    call theDumpBegins ( options )
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue

    if ( present(format) ) myFormat = format
    if ( index(myFormat, '*') > 0 ) &
      & myFormat = numNeedsFormat( one_c4*maxval(abs(array)), format )
    
    MyWidth = 3
    if ( present(width) ) MyWidth = width

    include 'dump3db.f9h'
  end subroutine DUMP_3D_DCOMPLEX

  ! ---------------------------------------------  DUMP_3D_DOUBLE  -----
  subroutine DUMP_3D_DOUBLE ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS, THESHAPE )
    double precision, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, J, K, L
    integer :: NumZeroRows
    double precision :: myFillValue
    character(len=64) :: myFormat
    integer :: nUnique
    integer :: MyWidth
    integer, dimension(MAXNUMELEMENTS) :: counts
    double precision, dimension(MAXNUMELEMENTS) :: elements
    myFormat = sdFormatDefault
    include 'dump3d.f9h'
    include 'dump3db.f9h'
  end subroutine DUMP_3D_DOUBLE

  ! --------------------------------------------  DUMP_3D_INTEGER  -----
  subroutine DUMP_3D_INTEGER ( ARRAY, NAME, &
    & FILLVALUE, FORMAT, WIDTH, LBOUND, OPTIONS, THESHAPE )
    integer, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many numbers per line (10)?
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, J, K, L
    integer :: NumZeroRows
    integer, dimension(3) :: which
    integer :: how_many

    character(len=64) :: MyFormat
    integer :: myFillValue
    integer :: nUnique
    integer :: MyWidth
    integer, dimension(MAXNUMELEMENTS) :: counts
    integer, dimension(MAXNUMELEMENTS) :: elements
    myFormat = 'places=' // INTPLACES ! To sneak places arg into call to output
    include 'dump3d.f9h'
    include 'dump3db.f9h'
  end subroutine DUMP_3D_INTEGER

  ! ---------------------------------------------  DUMP_3D_REAL  -----
  subroutine DUMP_3D_REAL ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS, THESHAPE )
    real, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    real, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: TheShape

    integer :: Base, I, J, K, L
    integer :: NumZeroRows
    real    :: myFillValue
    character(len=64) :: MyFormat
    integer :: myWidth
    integer :: nUnique
    integer, dimension(MAXNUMELEMENTS) :: counts
    real, dimension(MAXNUMELEMENTS) :: elements
    myFormat = sdFormatDefault
    include 'dump3d.f9h'
    include 'dump3db.f9h'
  end subroutine DUMP_3D_REAL

  ! ---------------------------------------------  DUMP_4D_DOUBLE  -----
  subroutine DUMP_4D_DOUBLE ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    double precision, intent(in) :: ARRAY(:,:,:,:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options

    integer :: J
    character(len=64) :: myFormat
    integer :: nUnique
    integer, dimension(MAXNUMELEMENTS) :: counts
    double precision, dimension(MAXNUMELEMENTS) :: elements
    myFormat = sdFormatDefault
    include 'dump4d.f9h'
    include 'dump4db.f9h'
  end subroutine DUMP_4D_DOUBLE

  ! ---------------------------------------------  DUMP_4D_REAL  -----
  subroutine DUMP_4D_REAL ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    real, intent(in) :: ARRAY(:,:,:,:)
    character(len=*), intent(in), optional :: NAME
    real, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options

    integer :: J
    character(len=64) :: myFormat
    integer :: nUnique
    integer, dimension(MAXNUMELEMENTS) :: counts
    real, dimension(MAXNUMELEMENTS) :: elements
    myFormat = sdFormatDefault
    include 'dump4d.f9h'
    include 'dump4db.f9h'
  end subroutine DUMP_4D_real

  ! -----------------------------------------------  DUMP_HASH_STR  -----
  subroutine DUMP_HASH_STR ( COUNTEMPTY, KEYS, VALUES, NAME, SEPARATOR, OPTIONS )
    ! Dumps a hash composed of two string lists: keys and values
    logical, intent(in)          :: COUNTEMPTY
    character(len=*), intent(in) :: KEYS
    character(len=*), intent(in) :: VALUES
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: SEPARATOR
    character(len=*), intent(in), optional :: options

    character( len=max(len(values), len(keys)) ) :: element
    integer :: J
    integer :: NumElements
    character(len=*), parameter :: COMMA = ','

    call theDumpBegins ( options )

    NumElements = NumStringElements(keys, countEmpty, &
     & separator)
    if ( NumElements == 0 ) then
      call empty ( name )
    else
      call name_and_size ( name, .false., NumElements )
      if ( present(name) .and. .not. mylaconic ) call newLine
      call output ( '(key)    =   (value)', advance='yes' )
      do j = 1, NumElements
        call GetStringElement( keys, element, j, countEmpty, separator)
        call output ( trim(element), advance='no' )
        call output( ' = ', advance='no' )
        call GetStringElement( values, element, j, countEmpty, separator)
        call output ( trim(element), advance='yes' )
      end do ! j
    end if
    call theDumpEnds
  end subroutine DUMP_HASH_STR

  ! -----------------------------------------------  DUMP_HASH_LOG  -----
  subroutine DUMP_HASH_LOG ( COUNTEMPTY, KEYS, VALUES, NAME, SEPARATOR, OPTIONS )
    ! Dumps a hash composed of two string lists: keys and values
    logical, intent(in)          :: COUNTEMPTY
    character(len=*), intent(in) :: KEYS
    logical, dimension(:), intent(in) :: VALUES
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: SEPARATOR
    character(len=*), intent(in), optional :: OPTIONS

    character( len=len(keys)) :: element
    integer :: J
    integer :: NumElements
    character(len=*), parameter :: COMMA = ','

    call theDumpBegins ( options )

    NumElements = NumStringElements(keys, countEmpty, &
     & separator)
    if ( NumElements == 0 ) then
      call empty ( name )
    else
      call name_and_size ( name, .false., NumElements )
      if ( present(name) .and. .not. mylaconic ) call newLine
      call output ( '(key)    =   (value)', advance='yes' )
      do j = 1, NumElements
        call GetStringElement( keys, element, j, countEmpty, separator)
        call output ( trim(element), advance='no' )
        call output( ' = ', advance='no' )
        call output ( values(j), advance='yes' )
      end do ! j
    end if
    call theDumpEnds
  end subroutine DUMP_HASH_LOG

  ! -----------------------------------------------  DUMP_STRLIST  -----
  subroutine DUMP_STRLIST ( STRING, NAME, FILLVALUE, OPTIONS, INSEPARATOR )
    ! Dumps a ','-separated string list, one item per lines
    ! (Unless it consists of multiple lines)
    character(len=*), intent(in) :: STRING
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: OPTIONS
    character(len=*), optional, intent(in) :: INSEPARATOR

    integer :: J
    integer :: NumElements
    character(len=len(string)) :: myFillValue
    character(len=1), parameter :: CR = ACHAR(13) ! Carriage return
    character(len=1), parameter :: LF = ACHAR(10) ! Line feed
    character(len=1) :: SEPARATOR
    logical, parameter :: COUNTEMPTY = .true.
    ! Executable
    if( index(string, CR) > 0 ) then
      call dump( (/ trim(string) /), name, fillvalue, options=options )
      return
    elseif( index(string, LF) > 0 ) then
      call dump( (/ trim(string) /), name, fillvalue, options=options )
      return
    endif

    call theDumpBegins ( options )
    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue

    SEPARATOR = ','
    if ( present(INSEPARATOR) ) SEPARATOR = INSEPARATOR
    
    NumElements = NumStringElements(string, countEmpty, &
     & separator)
    if ( NumElements == 0 ) then
      call empty ( name )
    else if ( NumElements == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( trim(string), advance='yes' )
    else
      call name_and_size ( name, myClean, NumElements )
      if ( present(name) .and. .not. mylaconic ) call newLine
      do j = 1, NumElements
        call GetStringElement(string, myFillValue, j, countEmpty, separator)
        call output ( trim(myFillValue), advance='yes' )
      end do ! j
    end if
    call theDumpEnds
  end subroutine DUMP_STRLIST

  ! ---------------------------------------------- DumpDates -----
  ! This family dumps dates (re)formatted according to dates_module
  subroutine dump_tai( taiDates, name, width, dateFormat, timeFormat )
    ! Dump an array of tai dates in whatever formats the user specifies
    ! Warning: does not bother with leap seconds
    ! By default we dump both date and time fields
    ! If dateFormat is present and blank or 'none', don't print date
    ! If timeFormat is present and blank or 'none', don't print time
    ! Args
    double precision, dimension(:)       :: taiDates ! tai93 (s)
    character(len=*), optional           :: name
    integer, intent(in), optional        :: WIDTH
    character(len=*), optional           :: dateFormat
    character(len=*), optional           :: timeFormat
    ! Internal variables
    character(len=16)                            :: date, time
    character(len=MAXUTCSTRLENGTH), dimension(size(taiDates)) &
      &                                          :: dates
    integer                                      :: error
    integer                                      :: i
    ! Executable
    if ( size(taidates) < 1 ) then
      call dump ( taiDates, name )
      return
    endif
    do i = 1, size(taiDates)
      dates(i) = tai93s2utc( taiDates(i) )
    enddo
    if ( present(dateFormat) .or. present(timeFormat) ) then
      do i = 1, size(taiDates)
        call splitDateTime( dates(i), error, date, time )
        if ( present(dateFormat) ) then
          if ( len_trim(dateFormat) < 1 .or. lowercase(dateFormat) == 'none' ) then
            date = ' '
          else
            date = ReformatDate( date, toForm=dateFormat )
          endif
        endif
        if ( present(timeFormat) ) then
          if ( len_trim(timeFormat) < 1 .or. lowercase(timeFormat) == 'none' ) then
            time = ' '
          else
            time = ReformatTime( time, Form=timeFormat )
          endif
        endif
        dates(i) = trim(date) // time
      enddo
    endif
    call dump( dates, name, width=width )
  end subroutine dump_tai
  
  ! ---------------------------------------------- DumpDumpOptions -----
  subroutine DumpDumpOptions( options )
    ! Show 
    ! (if current options supplied) actual dump, diff options
    ! (if arg is "?") available options and their meanings
    ! (if no arg) default options
    character(len=*), optional, intent(in) :: options
    character(len=1), parameter :: fillChar = '1' ! fill blanks with '. .'
     if ( .not. present(options) ) then
       call blanks(90, fillChar='-', advance='yes')
       call output(' ------------------------ Summary of automatic Dump, Diff options'      , advance='no')
       call output(' ------------------------ ', advance='yes')
       call outputNamedValue ( 'character printed between row, col id and data', aftersub, advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'default DIFF switches for CLEAN, TRIM, etc.', trim_safe(DEFAULTDIFFOPTIONS), advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'default DUMP switches for CLEAN, TRIM, etc.', trim_safe(DEFAULTDUMPOPTIONS), advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'print abs min, max, etc. when DIFF has RMS set TRUE?', DIFFRMSMEANSRMS, advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'skip dumping every element of a constant array?', DIFFRMSMEANSRMS, advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'what side to place headers when dumping tables', trim(DUMPTABLESIDE), advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'print stats all on one line?', STATSONONELINE, advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'pct output format', trim_safe(PCTFORMAT), advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'rms output format', trim_safe(RMSFORMAT), advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'numeric output format', trim_safe(SDFORMATDEFAULT), advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'complex output format', trim_safe(sdFormatDefaultCmplx), advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call blanks(90, fillChar='-', advance='yes')
     elseif ( index( options, "?" ) > 0 ) then
       call output( 'The meaning of options is determined by the presence or absence', advance='yes' )
       call output( 'of specific characters in the options string', advance='yes' )
       call output( 'If options is present and contains the following characters:', advance='yes' )
       call output( '(for dump or diff)', advance='yes' )
       call output( '  character         meaning', advance='yes' )
       call output( '     ---            -------', advance='yes' )
       call output( '      B              show Bandwidth, % of array that is non-zero', advance='yes' )
       call output( '      H              show rank, TheShape of array', advance='yes' )
       call output( '      L              laconic; skip printing name, size of array', advance='yes' )
       call output( '      N              show where NaNs and Infs are located', advance='yes' )
       call output( '      R              rms       -- min, max, etc.', advance='yes' )
       call output( '      b              table of % vs. amount of differences (pdf)', advance='yes' )
       call output( '      c              clean', advance='yes' )
       call output( '      g              gaps      ', advance='yes' )
       call output( '      l              collapse (last index)', advance='yes' )
       call output( '      r              ratios    -- min, max, etc. of difference ratios', advance='yes' )
       call output( '      s              stats     -- number of differences', advance='yes' )
       call output( '      p              transpose (for rank 2 arrays only)', advance='yes' )
       call output( '      t              trim      ', advance='yes' )
       call output( '      u              unique    ', advance='yes' )
       call output( '      w              wholearray', advance='yes' )
       call output( '      W[i]           wholearray, looping over ith index', advance='yes' )
       call output( '                     (for rank 3 and 4 arrays only)', advance='yes' )
       call output( '      1 or 2 or ..   ignored; calling routine is free to interpret', advance='yes' )
       call output( ' ', advance='yes' )
       call output( 'An exception is the behavior of w (wholearray):', advance='yes' )
       call output( 'if all {HNRblrs} are FALSE, i.e. unset, the whole array is dumped (or diffed)', advance='yes' )
       call output( 'if any are TRUE the whole array will be dumped only if', advance='yes' )
       call output( 'w is present or wholearray is set to TRUE', advance='yes' )
     else
       call outputNamedValue ( 'options', trim_safe(options), advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'thisIsADiff?', thisIsADiff, advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'myClean?', myClean, advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'myCollapse?', myCollapse, advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'myCyclic?', myCyclic, advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'myRMS?', myRMS, advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'myRatios?', myRatios, advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
       call outputNamedValue ( 'myWholeArray?', myWholeArray, advance='yes', &
         & fillChar=fillChar, before='* ', after=' *', tabn=4, tabc=62, taba=90 )
     endif
  end subroutine DumpDumpOptions

  ! -----------------------------------  dumpLists  -----
  ! Dump a 2d array as a set of lists
  ! Show names and related (numerical) values
  subroutine dumpLists_chars ( ARRAY, NAME, WIDTH, SEP, DELIMS )
    ! Args
    character(len=*), dimension(:,:), intent(in)   :: array
    character(len=*), intent(in), optional         :: NAME
    integer, optional, intent(in)                  :: width
    character(len=*), optional, intent(in)         :: sep
    character(len=*), optional, intent(in)         :: delims
    ! Internal variables
    character(len=8) :: intChars
    integer :: leftMargin
    integer :: line
    integer :: m
    integer :: myWidth ! how wide to make each list
    integer :: n
    integer :: n1
    integer :: n2
    integer :: nElemsPerLine
    integer :: nlines
    character(len=32) :: tabRange
    ! Executable
    call theDumpBegins ( options=' ' )
    if ( size(array, 2) < 1 ) then
      call empty(name)
    endif
    call name_and_size( name, clean=.false., size=size(array, 2) )
    m = size(array,1) ! num of elements in each list
    myWidth = (len(array(1,1))+2) * m + 3
    if ( present(width) ) myWidth = width
    leftMargin = max(4,ilog10(size(array,2)+1)) + 2
    nElemsPerLine = max( 1, (MAXLINELEN-leftMargin-1)/myWidth )
    nLines = 1 + (size(array,2)-1)/nElemsPerLine
    call writeIntsToChars( leftMargin, intChars )
    tabRange = intChars
    call writeIntsToChars( MAXLINELEN, intChars )
    tabRange = trim(tabRange) // '-' // intChars
    call writeIntsToChars( myWidth, intChars )
    tabRange = trim(tabRange) // '+' // intChars
    call setTabs( range=tabRange )
    n2 = 0
    do line=1, nLines
      n1 = n2 + 1
      n2 = min( n2 + nElemsPerLine, size(array,2) )
      call newLine
      call output ( n1, places=max(4,ilog10(size(array,2)+1)) , advance='no' )
      call output ( afterSub , advance='no' )
      do n=n1, n2
        call blanksToTab
        call outputList( array(:,n), sep, delims)
      enddo
    enddo
  end subroutine dumpLists_chars

  subroutine dumpLists_ints ( ARRAY, NAME, WIDTH, SEP, DELIMS )
    ! Args
    integer, dimension(:,:), intent(in)            :: array
    character(len=*), intent(in), optional         :: NAME
    character(len=*), optional, intent(in)         :: sep
    integer, optional, intent(in)                  :: width
    character(len=*), optional, intent(in)         :: delims
    ! Internal variables
    integer :: arrayMax
    character(len=8) :: intChars
    integer :: leftMargin
    integer :: line
    integer :: m
    integer :: myWidth ! how wide to make each list
    integer :: n
    integer :: n1
    integer :: n2
    integer :: nElemsPerLine
    integer :: nlines
    character(len=32) :: tabRange
    ! Executable
    call theDumpBegins ( options=' ' )
    if ( size(array, 2) < 1 ) then
      call empty(name)
    endif
    call name_and_size( name, clean=.false., size=size(array, 2) )
    m = size(array,1) ! num of elements in each list
    arrayMax = maxval(abs(array)) + 1
    myWidth = ( ilog10(arrayMax) + 1 ) * m + 3
    if ( present(width) ) myWidth = width
    leftMargin = max(4,ilog10(size(array,2)+1)) + 2
    nElemsPerLine = max( 1, (MAXLINELEN-leftMargin-1)/myWidth )
    nLines = 1 + (size(array,2)-1)/nElemsPerLine
    call writeIntsToChars( max(4,ilog10(size(array,2)+1)) + 4, intChars )
    tabRange = intChars
    call writeIntsToChars( MAXLINELEN, intChars )
    tabRange = trim(tabRange) // '-' // intChars
    call writeIntsToChars( myWidth, intChars )
    tabRange = trim(tabRange) // '+' // intChars
    call setTabs( range=tabRange )
    n2 = 0
    do line=1, nLines
      n1 = n2 + 1
      n2 = min( n2 + nElemsPerLine, size(array,2) )
      call newLine
      call output ( n1, places=max(4,ilog10(size(array,2)+1)) , advance='no' )
      call output ( afterSub , advance='no' )
      do n=n1, n2
        call blanksToTab
        call outputList( array(:,n), sep, delims)
      enddo
    enddo
    ! Reset tabs
    call newLine
    call resetTabs
  end subroutine dumpLists_ints

  ! -----------------------------------  dumpNamedValues  -----
  ! Another hash-like dump:
  ! Show names and related (numerical) values
  ! -----------------------------------  dumpNamedValues_DOUBLE  -----
  subroutine dumpNamedValues_DOUBLE ( VALUES, NAMES, FORMAT, WIDTH, OPTIONS )
    double precision, intent(in)                         :: values(:)
    character(len=*), intent(in), optional :: NAMES
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many pairs per line (1)?
    character(len=*), intent(in), optional :: options
    
    integer :: J, K, L
    integer :: MyWidth
    character(len=24) :: myName
    call theDumpBegins ( options )
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    if ( myClean ) call output(' ', advance='yes')
    l = 0
    do j=1, size(values), MyWidth
      do k=1, MyWidth
        call blanks(3, advance='no')
        l = l + 1
        if ( l <= size(values) ) then
          if ( present (names) ) then
            call GetStringElement(names, myName, l, .true.)
          else
            write(myName, *) 'double # ', l, ': '
          end if
          call output(myName,  advance='no')
          call blanks(3, advance='no')
          call output(values(l), format, advance='no')
        end if
      end do
      call output(' ', advance='yes')
    end do
    call theDumpEnds

  end subroutine dumpNamedValues_DOUBLE

  ! ----------------------------------  dumpNamedValues_INTEGER  -----
  subroutine dumpNamedValues_INTEGER ( VALUES, NAMES, FORMAT, WIDTH, OPTIONS )
    integer, intent(in)                         :: values(:)
    character(len=*), intent(in), optional :: NAMES
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many pairs per line (1)?
    character(len=*), intent(in), optional :: options
    
    integer :: J, K, L
    integer :: MyWidth
    character(len=24) :: myName
    call theDumpBegins ( options )
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    if ( myClean ) call output(' ', advance='yes')
    l = 0
    do j=1, size(values), MyWidth
      do k=1, MyWidth
        call blanks(3, advance='no')
        l = l + 1
        if ( l <= size(values) ) then
          if ( present (names) ) then
            call GetStringElement(names, myName, l, .true.)
          else
            write(myName, *) 'integer # ', l, ': '
          end if
          call output(myName,  advance='no')
          call blanks(3, advance='no')
          call output(values(l), format=format, advance='no')
        end if
      end do
      call output(' ', advance='yes')
    end do
    call theDumpEnds

  end subroutine dumpNamedValues_INTEGER

  ! -------------------------------------  dumpNamedValues_REAL  -----
  subroutine dumpNamedValues_REAL ( VALUES, NAMES, FORMAT, WIDTH, OPTIONS )
    real, intent(in)                         :: values(:)
    character(len=*), intent(in), optional :: NAMES
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many pairs per line (1)?
    character(len=*), intent(in), optional :: options
    
    integer :: J, K, L
    integer :: MyWidth
    character(len=24) :: myName
    call theDumpBegins ( options )
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    if ( myClean ) call output(' ', advance='yes')
    l = 0
    do j=1, size(values), MyWidth
      do k=1, MyWidth
        call blanks(3, advance='no')
        l = l + 1
        if ( l <= size(values) ) then
          if ( present (names) ) then
            call GetStringElement(names, myName, l, .true.)
          else
            write(myName, *) 'single # ', l, ': '
          end if
          call output(myName,  advance='no')
          call blanks(3, advance='no')
          call output(values(l), format, advance='no')
        end if
      end do
      call output(' ', advance='yes')
    end do
    call theDumpEnds

  end subroutine dumpNamedValues_REAL

 ! ---------------------------------------------  DUMPSUMS  -----
 ! This family of routines dumps the running sum:
 ! summed(i) == ( array(i) + summed(i-1) )
  subroutine DUMPSUMS_DOUBLE ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    double precision, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options
    ! Local variables
    integer                                  :: i
    double precision, dimension(size(array)) :: summed
    ! Executable
    if ( size(array) < 1 ) return
    summed(1) = array(1)
    if ( size(array) > 1 ) then
      do i=2, size(array)
        summed(i) = array(i) + summed(i-1)
      enddo
    endif
    call dump ( summed, NAME, &
      & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
  end subroutine DUMPSUMS_DOUBLE

  subroutine DUMPSUMS_INTEGER ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options
    ! Local variables
    integer                                  :: i
    integer, dimension(size(array)) :: summed
    ! Executable
    if ( size(array) < 1 ) return
    summed(1) = array(1)
    if ( size(array) > 1 ) then
      do i=2, size(array)
        summed(i) = array(i) + summed(i-1)
      enddo
    endif
    call dump ( summed, NAME, &
      & FILLVALUE, FORMAT, WIDTH, LBOUND, OPTIONS )
  end subroutine DUMPSUMS_INTEGER

  subroutine DUMPSUMS_REAL ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    real, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    real, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options
    ! Local variables
    integer                                  :: i
    real, dimension(size(array)) :: summed
    ! Executable
    if ( size(array) < 1 ) return
    summed(1) = array(1)
    if ( size(array) > 1 ) then
      do i=2, size(array)
        summed(i) = array(i) + summed(i-1)
      enddo
    endif
    call dump ( summed, NAME, &
      & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
  end subroutine DUMPSUMS_REAL

  ! -----------------------------------  dumpTable  -----
  ! This family of routines dumps a table of values and headers
  ! The headside determines on which of the 4 sides 
  ! ('top', 'left', 'right', 'bottom') we place the headers
  subroutine dumpTable_DOUBLE ( VALUES, HEADERS, HEADSIDE, FORMAT, FORMATS )
    double precision, dimension(:,:), intent(in) :: values
    character(len=*), dimension(:), intent(in)   :: HEADERS
    character(len=*), intent(in)                 :: HEADSIDE
    character(len=*), intent(in), optional       :: FORMAT
    character(len=*), dimension(:), intent(in), optional :: FORMATS
    include 'dumpTable.f9h'
  end subroutine dumpTable_DOUBLE

  subroutine dumpTable_INTEGER ( VALUES, HEADERS, HEADSIDE, FORMAT, FORMATS )
    integer, dimension(:,:), intent(in)          :: values
    character(len=*), dimension(:), intent(in)   :: HEADERS
    character(len=*), intent(in)                 :: HEADSIDE
    character(len=*), intent(in), optional       :: FORMAT
    character(len=*), dimension(:), intent(in), optional :: FORMATS
    include 'dumpTable.f9h'
  end subroutine dumpTable_INTEGER

  subroutine dumpTable_REAL ( VALUES, HEADERS, HEADSIDE, FORMAT, FORMATS )
    real, dimension(:,:), intent(in)             :: values
    character(len=*), dimension(:), intent(in)   :: HEADERS
    character(len=*), intent(in)                 :: HEADSIDE
    character(len=*), intent(in), optional       :: FORMAT
    character(len=*), dimension(:), intent(in), optional :: FORMATS
    include 'dumpTable.f9h'
  end subroutine dumpTable_REAL

 ! ---------------------------------------------  SELFDIFF_DOUBLE  -----
  subroutine SELFDIFF_DOUBLE ( ARRAY, NAME, &
    & FILLVALUE, WIDTH, FORMAT, waves, LBOUND, OPTIONS )
    ! dump the running increment == ( array(i) - array(i-1) )
    double precision, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WAVES
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options
    ! Local variables
    integer                                  :: i
    double precision, dimension(size(array)-1) :: increment
    integer, dimension(100) :: lengths
    logical :: myWaves
    integer :: nWaves
    ! Executable
    myWaves = .false.
    if ( present(waves) ) myWaves = waves
    if ( size(array) < 2 ) return
    do i=2, size(array)
      increment(i-1) = array(i) - array(i-1)
    enddo
    if ( myWaves ) then
      call halfWaves( increment, lengths, nWaves )
      if ( nWaves > 0 ) then
        call dump( lengths(1:nWaves), 'half-waves of ' // NAME )
      endif
      return
    endif
    call dump ( increment, NAME, &
      & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
  end subroutine SELFDIFF_DOUBLE

 ! ---------------------------------------------  SELFDIFF_INTEGER  -----
  subroutine SELFDIFF_INTEGER ( ARRAY, NAME, &
    & FILLVALUE, FORMAT, WIDTH, waves, LBOUND, OPTIONS )
    ! dump the running increment == ( array(i) - array(i-1) )
    integer, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WAVES
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options
    ! Local variables
    integer                                  :: i
    integer, dimension(size(array)-1) :: increment
    ! Executable
    if ( size(array) < 2 ) return
    do i=2, size(array)
      increment(i-1) = array(i) - array(i-1)
    enddo
    call dump ( increment, NAME, &
      & FILLVALUE, FORMAT, WIDTH, LBOUND, OPTIONS )
  end subroutine SELFDIFF_INTEGER

 ! ---------------------------------------------  SELFDIFF_REAL  -----
  subroutine SELFDIFF_REAL ( ARRAY, NAME, &
    & FILLVALUE, FORMAT, WIDTH, waves, LBOUND, OPTIONS )
    ! dump the running increment == ( array(i) - array(i-1) )
    real, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    real, intent(in), optional :: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WAVES
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options
    ! Local variables
    integer                                  :: i
    real, dimension(size(array)-1) :: increment
    ! logical, parameter :: unique = .false. 
    integer, dimension(100) :: lengths
    logical :: myWaves
    integer :: nWaves
    ! Executable
    myWaves = .false.
    if ( present(waves) ) myWaves = waves
    if ( size(array) < 2 ) return
    do i=2, size(array)
      increment(i-1) = array(i) - array(i-1)
    enddo
    if ( myWaves ) then
      call halfWaves( increment, lengths, nWaves )
      if ( nWaves > 0 ) then
        call dump( lengths(1:nWaves), 'half-waves of ' // NAME )
      endif
      return
    endif
    call dump ( increment, NAME, &
      & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
  end subroutine SELFDIFF_REAL

  ! --- Private procedures ---
  ! --- DumpCollapsedArray ---
  ! This family of subroutines dumps a lower-rank representation of
  ! an array
  subroutine DUMPCOLLAPSEDARRAY_1D_INTEGER (  array, name, fillValue, options )
    INTEGER, intent(in) :: ARRAY(:)
    character(len=*), intent(in) :: NAME
    INTEGER, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: options
    call output( 'Did not expect to dump collapsed 1d array', advance='yes' )
  end subroutine DUMPCOLLAPSEDARRAY_1D_INTEGER

  subroutine DUMPCOLLAPSEDARRAY_1D_DOUBLE (  array, name, fillValue, options )
    double precision, intent(in) :: ARRAY(:)
    character(len=*), intent(in) :: NAME
    double precision, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: options
    call output( 'Did not expect to dump collapsed 1d array', advance='yes' )
  end subroutine DUMPCOLLAPSEDARRAY_1D_DOUBLE

  subroutine DUMPCOLLAPSEDARRAY_1D_REAL (  array, name, fillValue, options )
    real, intent(in) :: ARRAY(:)
    character(len=*), intent(in) :: NAME
    real, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: options
    call output( 'Did not expect to dump collapsed 1d array', advance='yes' )
  end subroutine DUMPCOLLAPSEDARRAY_1D_REAL

  subroutine DUMPCOLLAPSEDARRAY_2D_DOUBLE (  array, name, fillValue, options )
    double precision, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in) :: NAME
    double precision, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: options
    ! For dumping lower-rank collapsed representations
    double precision, dimension(size(array, 1)) :: nums
    logical, dimension(size(array, 1))          :: logs

    ! dump numerical representation
    if ( index(COLLAPSEOPTIONS, 'num') > 0 ) then
      call collapse( array, nums, options=COLLAPSEOPTIONS )
      call dump( nums, name, fillvalue, options=options )
    endif
    if ( index(COLLAPSEOPTIONS, 'any') > 0 .or. &
      &  index(COLLAPSEOPTIONS, 'all') > 0 ) then
      call collapse( array, logs=logs, options=COLLAPSEOPTIONS )
      call dump( logs, name, options=options )
    endif
  end subroutine DUMPCOLLAPSEDARRAY_2D_DOUBLE

  subroutine DUMPCOLLAPSEDARRAY_2D_REAL (  array, name, fillValue, options )
    real, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in) :: NAME
    real, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: options
    ! For dumping lower-rank collapsed representations
    real, dimension(size(array, 1)) :: nums
    logical, dimension(size(array, 1))          :: logs

    ! dump numerical representation
    if ( index(COLLAPSEOPTIONS, 'num') > 0 ) then
      call collapse( array, nums, options=COLLAPSEOPTIONS )
      call dump( nums, name, fillvalue, options=options )
    endif
    if ( index(COLLAPSEOPTIONS, 'any') > 0 .or. &
      &  index(COLLAPSEOPTIONS, 'all') > 0 ) then
      call collapse( array, logs=logs, options=COLLAPSEOPTIONS )
      call dump( logs, name, options=options )
    endif
  end subroutine DUMPCOLLAPSEDARRAY_2D_REAL

  subroutine DUMPCOLLAPSEDARRAY_2D_INTEGER (  array, name, fillValue, options )
    INTEGER, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in) :: NAME
    INTEGER, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: options
    ! For dumping lower-rank collapsed representations
    INTEGER, dimension(size(array, 1)) :: nums
    logical, dimension(size(array, 1))          :: logs

    ! dump numerical representation
    if ( index(COLLAPSEOPTIONS, 'num') > 0 ) then
      call collapse( array, nums, options=COLLAPSEOPTIONS )
      call dump( nums, name, fillvalue, options=options )
    endif
    if ( index(COLLAPSEOPTIONS, 'any') > 0 .or. &
      &  index(COLLAPSEOPTIONS, 'all') > 0 ) then
      call collapse( array, logs=logs, options=COLLAPSEOPTIONS )
      call dump( logs, name, options=options )
    endif
  end subroutine DUMPCOLLAPSEDARRAY_2D_INTEGER

  subroutine DUMPCOLLAPSEDARRAY_3D_DOUBLE (  array, name, fillValue, options )
    double precision, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in) :: NAME
    double precision, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: options
    ! For dumping lower-rank collapsed representations
    double precision, dimension(size(array, 1),size(array, 2)) :: nums
    logical, dimension(size(array, 1),size(array, 2))          :: logs

    ! dump numerical representation
    if ( index(COLLAPSEOPTIONS, 'num') > 0 ) then
      call collapse( array, nums, options=COLLAPSEOPTIONS )
      call dump( nums, name, fillvalue, options=options )
    endif
    if ( index(COLLAPSEOPTIONS, 'any') > 0 .or. &
      &  index(COLLAPSEOPTIONS, 'all') > 0 ) then
      call collapse( array, logs=logs, options=COLLAPSEOPTIONS )
      call dump( logs, name, options=options )
    endif
  end subroutine DUMPCOLLAPSEDARRAY_3D_DOUBLE

  subroutine DUMPCOLLAPSEDARRAY_3D_REAL (  array, name, fillValue, options )
    real, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in) :: NAME
    real, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: options
    ! For dumping lower-rank collapsed representations
    real, dimension(size(array, 1),size(array, 2)) :: nums
    logical, dimension(size(array, 1),size(array, 2))          :: logs

    ! dump numerical representation
    if ( index(COLLAPSEOPTIONS, 'num') > 0 ) then
      call collapse( array, nums, options=COLLAPSEOPTIONS )
      call dump( nums, name, fillvalue, options=options )
    endif
    if ( index(COLLAPSEOPTIONS, 'any') > 0 .or. &
      &  index(COLLAPSEOPTIONS, 'all') > 0 ) then
      call collapse( array, logs=logs, options=COLLAPSEOPTIONS )
      call dump( logs, name, options=options )
    endif
  end subroutine DUMPCOLLAPSEDARRAY_3D_REAL

  subroutine DUMPCOLLAPSEDARRAY_3D_INTEGER (  array, name, fillValue, options )
    INTEGER, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in) :: NAME
    INTEGER, intent(in), optional :: FILLVALUE
    character(len=*), intent(in), optional :: options
    ! For dumping lower-rank collapsed representations
    INTEGER, dimension(size(array, 1),size(array, 2)) :: nums
    logical, dimension(size(array, 1),size(array, 2))          :: logs

    ! dump numerical representation
    if ( index(COLLAPSEOPTIONS, 'num') > 0 ) then
      call collapse( array, nums, options=COLLAPSEOPTIONS )
      call dump( nums, name, fillvalue, options=options )
    endif
    if ( index(COLLAPSEOPTIONS, 'any') > 0 .or. &
      &  index(COLLAPSEOPTIONS, 'all') > 0 ) then
      call collapse( array, logs=logs, options=COLLAPSEOPTIONS )
      call dump( logs, name, options=options )
    endif
  end subroutine DUMPCOLLAPSEDARRAY_3D_INTEGER

  !
  subroutine FILTEREDDIFF_1D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    double precision, intent(in) :: ARRAY1(:)
    character(len=*), intent(in) :: NAME1
    double precision, intent(in) :: ARRAY2(:)
    character(len=*), intent(in) :: NAME2
    double precision, intent(in):: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    double precision, dimension(size(array1)) :: filtered1
    double precision, dimension(size(array2)) :: filtered2
    double precision :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FILTEREDDIFF_1D_DOUBLE

  subroutine FILTEREDDIFF_1D_INTEGER ( IARRAY1, NAME1, IARRAY2, NAME2, &
    & IFILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, intent(in) :: IARRAY1(:)
    character(len=*), intent(in) :: NAME1
    integer, intent(in) :: IARRAY2(:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: IFILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    real, dimension(size(iarray1)) :: array1
    real, dimension(size(iarray2)) :: array2
    real :: fillValue
    ! So we don't have to write an integer-version of allstats
    array1 = iarray1
    array2 = iarray2
    fillValue = iFillValue
    call DIFF ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
  end subroutine FILTEREDDIFF_1D_INTEGER

  subroutine FILTEREDDIFF_1D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    real, intent(in) :: ARRAY1(:)
    character(len=*), intent(in) :: NAME1
    real, intent(in) :: ARRAY2(:)
    character(len=*), intent(in) :: NAME2
    real, intent(in):: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    real, dimension(size(array1)) :: filtered1
    real, dimension(size(array2)) :: filtered2
    real :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FILTEREDDIFF_1D_REAL

  subroutine FILTEREDDIFF_2D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    double precision, intent(in) :: ARRAY1(:,:)
    character(len=*), intent(in) :: NAME1
    double precision, intent(in) :: ARRAY2(:,:)
    character(len=*), intent(in) :: NAME2
    double precision, intent(in):: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options
    !
    double precision, dimension(product(shape(array1))) :: filtered1
    double precision, dimension(product(shape(array2))) :: filtered2
    double precision :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FILTEREDDIFF_2D_DOUBLE

  subroutine FILTEREDDIFF_2D_INTEGER ( IARRAY1, NAME1, IARRAY2, NAME2, &
    & IFILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, intent(in) :: IARRAY1(:,:)
    character(len=*), intent(in) :: NAME1
    integer, intent(in) :: IARRAY2(:,:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: IFILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    real, dimension(size(iarray1,1), size(iarray1,2)) :: array1
    real, dimension(size(iarray1,1), size(iarray1,2)) :: array2
    real :: fillValue
    ! So we don't have to write an integer-version of allstats
    array1 = iarray1
    array2 = iarray2
    fillValue = iFillValue
    call DIFF ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
  end subroutine FILTEREDDIFF_2D_INTEGER

  subroutine FILTEREDDIFF_2D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    real, intent(in) :: ARRAY1(:,:)
    character(len=*), intent(in) :: NAME1
    real, intent(in) :: ARRAY2(:,:)
    character(len=*), intent(in) :: NAME2
    real, intent(in):: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options
    !
    real, dimension(product(shape(array1))) :: filtered1
    real, dimension(product(shape(array2))) :: filtered2
    real :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FILTEREDDIFF_2D_REAL

  subroutine FILTEREDDIFF_3D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    double precision, intent(in) :: ARRAY1(:,:,:)
    character(len=*), intent(in) :: NAME1
    double precision, intent(in) :: ARRAY2(:,:,:)
    character(len=*), intent(in) :: NAME2
    double precision, intent(in):: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options

    double precision, dimension(product(shape(array1))) :: filtered1
    double precision, dimension(product(shape(array2))) :: filtered2
    double precision :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FILTEREDDIFF_3D_DOUBLE

  subroutine FILTEREDDIFF_3D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    real, intent(in) :: ARRAY1(:,:,:)
    character(len=*), intent(in) :: NAME1
    real, intent(in) :: ARRAY2(:,:,:)
    character(len=*), intent(in) :: NAME2
    real, intent(in):: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options

    real, dimension(product(shape(array1))) :: filtered1
    real, dimension(product(shape(array2))) :: filtered2
    real :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FILTEREDDIFF_3D_REAL

  subroutine FILTEREDDIFF_4D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    double precision, intent(in) :: ARRAY1(:,:,:,:)
    character(len=*), intent(in) :: NAME1
    double precision, intent(in) :: ARRAY2(:,:,:,:)
    character(len=*), intent(in) :: NAME2
    double precision, intent(in):: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options

    double precision, dimension(product(shape(array1))) :: filtered1
    double precision, dimension(product(shape(array2))) :: filtered2
    double precision :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FILTEREDDIFF_4D_DOUBLE

  subroutine FILTEREDDIFF_4D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, WIDTH, FORMAT, LBOUND, OPTIONS )
    real, intent(in) :: ARRAY1(:,:,:,:)
    character(len=*), intent(in) :: NAME1
    real, intent(in) :: ARRAY2(:,:,:,:)
    character(len=*), intent(in) :: NAME2
    real, intent(in):: FILLVALUE
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options

    real, dimension(product(shape(array1))) :: filtered1
    real, dimension(product(shape(array2))) :: filtered2
    real :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FILTEREDDIFF_4D_REAL

  ! ---------------------------------------------  restoreDumpConfig  -----
  ! Restore default values for dump settings
  subroutine RESTOREDUMPCONFIG
    DEFAULTDIFFOPTIONS        = ' '
    DEFAULTDUMPOPTIONS        = ' '
    DEFAULTMAXLON             = 128
    DUMPTABLESIDE             = 'top'
    DIFFRMSMEANSRMS           = .false.
    DONTDUMPIFALLEQUAL        = .true.
    FILTERFILLSFROMRMS        = .false.
    PRINTFILLVALUE            = .true.
    PRINTNAMEIFDIFF           = .true.
    STATSONONELINE            = .true.
    COLLAPSEOPTIONS           = 'num[+]all[+]'
    INTPLACES                 = '6' ! how many places
    PCTFORMAT                 = '*' ! * means default format
    RMSFORMAT                 = '*' ! * means default format
    SDFORMATDEFAULT           = '(1pg14.6)'
  end subroutine RESTOREDUMPCONFIG

  ! ---------------------------------------------  arrayShapeToString  -----
  function arrayShapeToString ( arrayShape ) result ( string )
    ! Given an array of integers return the shape as a string
    ! E.g., given (/4,2,6/) return '4*2*6'
    integer, dimension(:), intent(in) :: arrayShape
    character(len=16) :: string
    ! Internal variables
    integer :: i
    ! Executable
    string = ' '
    if ( size(arrayShape) < 1 ) return
    do i=1, size( arrayshape )
      string = catLists( string, arrayShape(i), inseparator='*' )
    enddo
    
  end function arrayShapeToString
  ! ------------------------------------------------------  Empty  -----
  subroutine Empty ( Name )
    character(len=*), intent(in), optional :: Name

    if ( present(name) ) then
      call output ( name , advance='no' )
      call output ( ' is ' , advance='no' )
      nameHasBeenPrinted = .true.
    end if
    call output ( 'empty', advance='yes' )

  end subroutine Empty

  ! -----------------------------------------------------  ILOG10  -----
  integer function ILOG10(int)
    integer, intent(in) :: int
    ilog10=nint(log10(real(int)))
  end function ILOG10

  ! ----------------------------------------------  Name_And_Size  -----
  subroutine Name_And_Size ( NAME, CLEAN, SIZE, THESHAPE )
    character(len=*), intent(in), optional :: Name
    logical, intent(in) :: Clean
    integer, intent(in) :: Size
    character(len=*), intent(in), optional :: TheShape

    if ( present(name) .and. .not. myLaconic ) then
      if ( len_trim(name) < 1 ) return
      if ( .not. nameHasBeenPrinted ) then
        call output ( name , advance='no' )
        if ( present(theShape) ) call output ( theShape , advance='no' )
      end if
      if ( clean ) then 
        call output ( trim(" \ ") ) ! This goofiness is to outwit an incorrect
                                    ! Intel compiler.
        call output ( size , advance='no' )
      end if
      if ( size == 1 ) call output ( ' ' , advance='no' )
      nameHasBeenPrinted = .true.
    end if

  end subroutine Name_And_Size
  
  ! -------------------- PrintIt ---------------------
  ! This family of subroutines exists only so that we can generically call
  ! output with either a numeric arg or a character string, 
  ! trimming if the latter
  subroutine PrintIt_char ( it, format )
    character(len=*) :: it
    character(len=*), intent(in), optional :: FORMAT
    call output ( trim(it), advance='no' )
  end subroutine PrintIt_char

  subroutine PrintIt_int ( it, format )
    integer :: it
    character(len=*), intent(in), optional :: FORMAT
    call output ( it, format=format, advance='no' )
  end subroutine PrintIt_int

  subroutine PrintIt_real ( it, format )
    real :: it
    character(len=*), intent(in), optional :: FORMAT
    call output ( it, format=format, advance='no' )
  end subroutine PrintIt_real

  subroutine PrintIt_double ( it, format )
    double precision :: it
    character(len=*), intent(in), optional :: FORMAT
    call output ( it, format=format, advance='no' )
  end subroutine PrintIt_double

  subroutine PrintIt_complex ( it, format )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk) :: it
    character(len=*), intent(in), optional :: FORMAT
    call output ( it, format=format, advance='no' )
  end subroutine PrintIt_complex

  subroutine PrintIt_dcomplex ( it, format )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk) :: it
    character(len=*), intent(in), optional :: FORMAT
    call output ( it, format=format, advance='no' )
  end subroutine PrintIt_dcomplex

  ! ----------------------------------------------  printPercentages  -----
  ! Prints a nicely-formatted summary of equal, unequal, etc.
  ! using output
  subroutine printPercentages ( name, equal, unequal )
    character(len=*), intent(in), optional :: Name
    integer, intent(in) :: equal
    integer, intent(in) :: unequal
    if ( equal+unequal < 1 ) return
    myPCTFormat  = DEFAULTPCTFORMAT
    if ( PCTFORMAT /= '*' ) myPCTFormat = PCTFORMAT
    if ( present(name) ) call output ( trim(name), advance='no' )
    call blanks( 1, advance='no' )
    call output( fillvaluerelation, advance='no' )
    call output ( ', !', advance='no' )
    call output( fillvaluerelation, advance='no' )
    call output ( ' (%) ', advance='no' )
    call output ( equal, advance='no' )
    call output ( ': ', advance='no' )
    call output ( unequal, advance='no' )
    if ( .not. STATSONONELINE ) then
      call newline
      call blanks(10)
    endif
    call output ( '( ', advance='no' )
    call output ( 100*equal/(equal+unequal+0.), format = myPCTFormat, advance='no' )
    call output ( ': ', advance='no' )
    call output ( 100*unequal/(equal+unequal+0.), format = myPCTFormat, advance='no' )
    call output ( ' )', advance='yes' )
  end subroutine printPercentages

  ! ----------------------------------------------  printRMSetc  -----
  ! This family of routines prints a nicely-formatted list of min, max, etc.
  ! using output
  subroutine printRMSetc_double ( Name, min, max, rms, mean  )
    character(len=*), intent(in), optional :: Name
    double precision, intent(in) :: min
    double precision, intent(in) :: max
    double precision, intent(in) :: rms
    double precision, intent(in), optional :: mean
    !
    character(len=16) :: originalSDFormat
    !
    originalSDFormat = outputOptions%sdFormatDefault
    outputOptions%sdFormatDefault = rmsFormat
    include 'printRMSetc.f9h'
    outputOptions%sdFormatDefault = originalSDFormat
  end subroutine printRMSetc_double

  subroutine printRMSetc_real ( Name, min, max, rms, mean  )
    character(len=*), intent(in), optional :: Name
    real, intent(in) :: min
    real, intent(in) :: max
    real, intent(in) :: rms
    real, intent(in), optional :: mean
    character(len=16) :: originalSDFormat
    !
    originalSDFormat = outputOptions%sdFormatDefault
    outputOptions%sdFormatDefault = rmsFormat
    include 'printRMSetc.f9h'
    outputOptions%sdFormatDefault = originalSDFormat
  end subroutine printRMSetc_real

  subroutine printRMSetc_int ( Name, min, max, rms, mean  )
    character(len=*), intent(in), optional :: Name
    integer, intent(in) :: min
    integer, intent(in) :: max
    real, intent(in) :: rms
    real, intent(in), optional :: mean
    include 'printRMSetc.f9h'
  end subroutine printRMSetc_int

  ! This family of subroutines print subscripts to the left
  ! of each dumped row, sometimes noting that repeated lines
  ! that have been omitted for brevity
  ! ----------------------------------------------  Say_Fill_Char  -----
  subroutine Say_Fill_Char ( Subs, NumZeroRows, Fill, Inc, Format  )
    character(len=*), intent(in) :: Fill
    include 'Say_Fill.f9h'
  end subroutine Say_Fill_Char

  ! --------------------------------------------  Say_Fill_Complex  -----
  subroutine Say_Fill_Complex ( Subs, NumZeroRows, Fill, Inc, Format )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: Fill
    include 'Say_Fill.f9h'
  end subroutine Say_Fill_Complex

  ! --------------------------------------------  Say_Fill_Dcomplex  -----
  subroutine Say_Fill_Dcomplex ( Subs, NumZeroRows, Fill, Inc, Format )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: Fill
    include 'Say_Fill.f9h'
  end subroutine Say_Fill_Dcomplex

  ! --------------------------------------------  Say_Fill_Double  -----
  subroutine Say_Fill_Double ( Subs, NumZeroRows, Fill, Inc, Format )
    double precision, intent(in) :: Fill
    include 'Say_Fill.f9h'
  end subroutine Say_Fill_Double

  ! -----------------------------------------------  Say_Fill_Int  -----
  subroutine Say_Fill_Int ( Subs, NumZeroRows, Fill, Inc, Format )
    integer, intent(in) :: Fill
    include 'Say_Fill.f9h'
  end subroutine Say_Fill_Int

  ! ----------------------------------------------  Say_Fill_Real  -----
  subroutine Say_Fill_Real ( Subs, NumZeroRows, Fill, Inc, Format )
    real, intent(in) :: Fill
    include 'Say_Fill.f9h'
  end subroutine Say_Fill_Real

  ! ----------------------------------------------------  Say_Subs -----
  subroutine Say_Subs ( Subs, NumZeroRows )
    integer, intent(in) :: Subs(:)
    integer, intent(in) :: NumZeroRows
    call say_subs_only ( subs )
    call output ( ' ' , advance='no' )
    call output ( numZeroRows , advance='no' )
    call output ( ' lines of ', advance='no' )
  end subroutine Say_Subs

  ! -----------------------------------------------  Say_Subs_Only -----
  subroutine Say_Subs_Only ( Subs )
    integer, intent(in) :: Subs(:)
    integer :: I
    do i = 1, size(subs), 2
      call output ( subs(i), places=max(4,ilog10(subs(i+1))+1) , advance='no' )
    end do
    call output ( afterSub , advance='no' )
  end subroutine Say_Subs_Only
  
  subroutine showColumnNums( lineLength, skip )
    ! show column numbers
    ! usu. printed below other data to help guide the eye as to where
    ! stuff actually is
    ! E.g.,
    !   TTTT  TT  T
    !       FF  FF FFFF
    !   12345678901234567890          <-- These are the 
    !            1         2          <-- two lines we print here
    ! Args:
    integer, intent(in)           :: lineLength ! How many per line
    integer, intent(in), optional :: skip       ! start numbering after skip
    ! Internal variables
    integer                          :: bloc
    character(len=10), dimension(2)  :: line
    integer                          :: multiple
    integer                          :: mySkip
    integer                          :: numBlocs
    ! Executable
    numBlocs = 1 + (lineLength-1)/10
    mySkip = 0
    if ( present(skip) ) mySkip = skip
    line(1) = '1234567890'
!   line(2) = '         1'
    do multiple=1, 2
      if ( mySkip > 0 ) call blanks(mySkip)
      do bloc = 1, numBlocs
        write(line(2), '(i10)' ) bloc
        call output( line(multiple), advance='no' )
      enddo
      call newline
    enddo
  end subroutine showColumnNums
  
  function snipOption( options, particular ) result ( snipped )
    ! Snip from options a particular option commanding dump or diff
    ! Args:
    character(len=*), optional, intent(in) :: options
    character(len=1), optional, intent(in) :: particular
    character(len=16)                      :: snipped
    ! Executable
    snipped = ' '
    if ( .not. present(options) .or. .not. present(particular) ) return
    snipped = delete( options, particular )
  end function snipOption
  
  subroutine theDumpBegins(options)
    character(len=*), intent(in), optional :: options
    nameHasBeenPrinted = .false.
    stampOptions%neverStamp = .true. ! So we don't interrupt tables of numbers
    myBandwidth  = theDefault('bandwidth') ! .false.
    myClean      = theDefault('clean') ! .false.
    myCollapse   = theDefault('collapse') ! .false.
    myCyclic     = theDefault('cyclic') ! .false.
    myDirect     = theDefault('direct') ! .false.
    myGaps       = theDefault('gaps')
    myLaconic    = theDefault('laconic')
    myNaNs       = theDefault('nans')
    myRMS        = theDefault('rms')   ! .false.
    MyRatios     = theDefault('ratios')   ! .false.
    myShape      = theDefault('shape')  ! .false.
    myStats      = theDefault('stat')  ! .false.
    myTable      = theDefault('table')  ! .false.
    myTranspose  = theDefault('transpose')  ! .false.
    myTrim       = theDefault('trim')  ! .false.
    myUnique     = theDefault('unique')
    myWholeArray = theDefault('wholearray')
    if ( present(options) ) then
      myBandwidth   =   index( options, dopt_bandwidth  ) > 0
      myClean       =   index( options, dopt_clean      ) > 0
      myCollapse    =   index( options, dopt_collapse   ) > 0
      myCyclic      =   index( options, dopt_cyclic     ) > 0
      myGaps        =   index( options, dopt_gaps       ) > 0
      myLaconic     =   index( options, dopt_laconic    ) > 0
      myNaNs        =   index( options, dopt_NaNs       ) > 0
      myRatios      =   index( options, dopt_ratios     ) > 0
      myRMS         =   index( options, dopt_rms        ) > 0
      myShape       =   index( options, dopt_shape      ) > 0
      myStats       =   index( options, dopt_stats      ) > 0
      myTable       =   index( options, dopt_table      ) > 0
      myTranspose   =   index( options, dopt_transpose  ) > 0
      myTrim        =   index( options, dopt_trim       ) > 0
      myUnique      =   index( options, dopt_unique     ) > 0
      myWholeArray  = ( index( options, dopt_wholearray ) > 0 )
    endif
    myWholeArray = myWholeArray .or. &
      & .not. (myBandwidth.or. myCollapse .or. myRatios .or. myRMS .or. myShape .or. myStats &
      & .or. myTable .or. myNaNs)
    onlyWholeArray = myWholeArray .and. &
      & .not. (myBandwidth.or. myCollapse .or. myRatios .or. myRMS .or. myShape .or. myStats &
      & .or. myTable .or. myNaNs)
    nameHasBeenPrinted = nameHasBeenPrinted .or. myLaconic
  end subroutine theDumpBegins

  subroutine theDumpEnds
    stampOptions%neverStamp = .false.
  end subroutine theDumpEnds

  function theDefault( code ) result ( isit )
    ! Return the default value for a given code, e.g. 'clean'
    ! Args
    character(len=*), intent(in)    :: code
    logical                         :: isit
    !
    character(len=8) :: defaultString
    ! Executable
    if ( thisIsADiff ) then
      defaultString = DEFAULTDIFFOPTIONS
    else
      defaultString = DEFAULTDUMPOPTIONS
    endif
    select case ( code )
    case ('clean')
      isit = index( defaultstring, dopt_clean      ) > 0
    case ('cyclic')
      isit = index( defaultstring, dopt_cyclic     ) > 0
    case ('direct')
      isit = index( defaultstring, dopt_collapse   ) > 0
    case ('gaps')
      isit = index( defaultstring, dopt_gaps       ) > 0
    case ('laconic')
      isit = index( defaultstring, dopt_laconic    ) > 0
    case ('nans')
      isit = index( defaultstring, dopt_NaNs       ) > 0
    case ('rms')
      isit = index( defaultstring, dopt_rms        ) > 0
    case ('shape')
      isit = index( defaultstring, dopt_shape      ) > 0
    case ('stat')
      isit = index( defaultstring, dopt_stats      ) > 0
    case ('table')
      isit = index( defaultstring, dopt_table      ) > 0
    case ('transpose')
      isit = index( defaultstring, dopt_transpose  ) > 0
    case ('trim')
      isit = index( defaultstring, dopt_trim       ) > 0
    case ('unique')
      isit = index( defaultstring, dopt_unique     ) > 0
    case ('wholearray')
      isit = index( defaultstring, dopt_wholearray ) > 0

    case default
      isit = .false.
    end select
  end function theDefault
  
  subroutine UNFILTEREDDIFF_1D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: ARRAY1(:)
    character(len=*), intent(in) :: NAME1
    real(rk), intent(in) :: ARRAY2(:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UNFILTEREDDIFF_1D_DOUBLE

  subroutine UNFILTEREDDIFF_1D_INTEGER ( IARRAY1, NAME1, IARRAY2, NAME2, &
    & WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, intent(in) :: IARRAY1(:)
    character(len=*), intent(in) :: NAME1
    integer, intent(in) :: IARRAY2(:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    real, dimension(size(iarray1)) :: array1
    real, dimension(size(iarray2)) :: array2
    ! So we don't have to write an integer-version of allstats
    array1 = iarray1
    array2 = iarray2
    call DIFF ( ARRAY1, NAME1, ARRAY2, NAME2, &
      & WIDTH=WIDTH, FORMAT=FORMAT, &
      & LBOUND=LBOUND, OPTIONS=OPTIONS )
  end subroutine UNFILTEREDDIFF_1D_INTEGER

  subroutine UNFILTEREDDIFF_1D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: ARRAY1(:)
    character(len=*), intent(in) :: NAME1
    real(rk), intent(in) :: ARRAY2(:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UNFILTEREDDIFF_1D_REAL

  subroutine UNFILTEREDDIFF_2D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: ARRAY1(:,:)
    character(len=*), intent(in) :: NAME1
    real(rk), intent(in) :: ARRAY2(:,:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UNFILTEREDDIFF_2D_DOUBLE

  subroutine UNFILTEREDDIFF_2D_INTEGER ( IARRAY1, NAME1, IARRAY2, NAME2, &
    & WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, intent(in) :: IARRAY1(:,:)
    character(len=*), intent(in) :: NAME1
    integer, intent(in) :: IARRAY2(:,:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    character(len=*), intent(in), optional :: options

    real, dimension(size(iarray1,1), size(iarray1,2)) :: array1
    real, dimension(size(iarray1,1), size(iarray1,2)) :: array2
    ! So we don't have to write an integer-version of allstats
    array1 = iarray1
    array2 = iarray2
    call DIFF ( ARRAY1, NAME1, ARRAY2, NAME2, &
      & WIDTH=WIDTH, FORMAT=FORMAT, &
      & LBOUND=LBOUND, OPTIONS=OPTIONS )
  end subroutine UNFILTEREDDIFF_2D_INTEGER

  subroutine UNFILTEREDDIFF_2D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: ARRAY1(:,:)
    character(len=*), intent(in) :: NAME1
    real(rk), intent(in) :: ARRAY2(:,:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UNFILTEREDDIFF_2D_REAL

  subroutine UNFILTEREDDIFF_3D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: ARRAY1(:,:,:)
    character(len=*), intent(in) :: NAME1
    real(rk), intent(in) :: ARRAY2(:,:,:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UNFILTEREDDIFF_3D_DOUBLE

  subroutine UNFILTEREDDIFF_3D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: ARRAY1(:,:,:)
    character(len=*), intent(in) :: NAME1
    real(rk), intent(in) :: ARRAY2(:,:,:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UNFILTEREDDIFF_3D_REAL

  subroutine UNFILTEREDDIFF_4D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: ARRAY1(:,:,:,:)
    character(len=*), intent(in) :: NAME1
    real(rk), intent(in) :: ARRAY2(:,:,:,:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UNFILTEREDDIFF_4D_DOUBLE

  subroutine UNFILTEREDDIFF_4D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & WIDTH, FORMAT, LBOUND, OPTIONS )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: ARRAY1(:,:,:,:)
    character(len=*), intent(in) :: NAME1
    real(rk), intent(in) :: ARRAY2(:,:,:,:)
    character(len=*), intent(in) :: NAME2
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: LBOUND
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UNFILTEREDDIFF_4D_REAL

  logical function uniqueonly ( options )
    character(len=*), intent(in), optional :: options
    ! Executable
    uniqueonly = .true.
    if ( .not. present(options) ) return
    uniqueonly = uniqueonly .and. &
      & .not. any( indexes(options, (/'w','s','r'/)) > 0 )

  end function uniqueonly

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module DUMP_0

! $Log$
! Revision 2.131  2014/08/06 23:02:21  vsnyder
! Use kind C_Int16_t from ISO_C_Bindings instead of INTEGER*2
!
! Revision 2.130  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.129  2013/09/26 15:25:41  pwagner
! Added 2-byte integer dumps
!
! Revision 2.128  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.127  2013/06/28 18:08:38  pwagner
! Note if all logical elements equal
!
! Revision 2.126  2013/06/19 23:14:39  pwagner
! Remove more unused stuff; give default value to myDirect
!
! Revision 2.125  2013/01/10 00:18:43  pwagner
! Dumps with -N just show where NaNs are
!
! Revision 2.124  2012/09/11 21:10:24  pwagner
! Requires 'N' option to show where NaNs, Infs are located
!
! Revision 2.123  2012/07/05 23:47:31  pwagner
! Added restoreDumpConfig
!
! Revision 2.122  2012/06/22 20:25:52  pwagner
! Specify advance arg because we may now set default to 'yes'
!
! Revision 2.121  2012/06/13 23:59:37  pwagner
! dumpDumpOptions optionally dumps available diff/dump options
!
! Revision 2.120  2012/01/09 22:25:55  pwagner
! Distinguish 'r' option to print rms of ratios and 'R' option for rms of values
!
! Revision 2.119  2011/12/15 01:47:45  pwagner
! Accepts W[i] option; deletes more prolix notices of rank reduction
!
! Revision 2.118  2011/12/07 01:14:44  pwagner
! Added option to show bandwidth of banded arrays
!
! Revision 2.117  2011/11/11 00:30:22  vsnyder
! Simplify notice of rank reduction
!
! Revision 2.116  2011/07/26 20:40:24  pwagner
! Added 4d diffs, too
!
! Revision 2.115  2011/07/23 00:16:25  vsnyder
! Use (1.,0.) instead of CMPLX(1.,0.) to avoid Intel gripe
!
! Revision 2.114  2011/07/15 23:23:43  pwagner
! Can now dump 4d arrays, with certain restrictions
!
! Revision 2.113  2011/07/12 00:15:01  pwagner
! Improved dumps; format option now more flexible; complex arrays parallel real ones
!
! Revision 2.112  2011/06/23 17:27:41  pwagner
! Added function to difference args with option to supply fillvalue or period
!
! Revision 2.111  2011/04/26 20:54:12  pwagner
! Can now dump dates
!
! Revision 2.110  2011/04/20 22:27:09  pwagner
! Fixed long-standing bug in dumping 2d char array
!
! Revision 2.109  2011/04/18 19:12:54  vsnyder
! Turn on nameHasBeenPrinted in more places
!
! Revision 2.108  2011/04/04 23:08:27  pwagner
! May diff 2d integer arrays
!
! Revision 2.107  2011/03/08 00:04:40  pwagner
! Skip printing row of zeros only if multiple
!
! Revision 2.106  2011/02/25 18:50:26  pwagner
! Added optional width arg to char array dumps
!
! Revision 2.105  2011/02/05 01:01:41  pwagner
! Added new dopt_ dump options; transpose and collapse work better
!
! Revision 2.104  2011/01/20 01:16:01  pwagner
! New '-l' option dumps collapsed representation of higher ranked array
!
! Revision 2.103  2011/01/04 00:48:26  pwagner
! DEFAULTMAXLON can now set max width of 1d char array dumps
!
! Revision 2.102  2010/10/14 18:44:01  pwagner
! Can now dump lists
!
! Revision 2.101  2010/02/04 23:05:39  vsnyder
! Remove USE and declaration for unreferenced names.
! Declare a kind parameter in unfiltereddiff*
!
! Revision 2.100  2010/01/29 21:08:21  pwagner
! gave myFillValue a value to stop lf95 complaints
!
! Revision 2.99  2009/11/20 01:15:11  pwagner
! Decided against using allFinite
!
! Revision 2.98  2009/11/20 01:12:50  pwagner
! Added new option H or 'shape' to just show array rank, Shape
!
! Revision 2.97  2009/10/30 23:02:41  pwagner
! Should not double-print name if only whole array diff
!
! Revision 2.96  2009/10/26 18:53:33  pwagner
! Fixed bug only NAG found
!
! Revision 2.95  2009/10/19 17:33:26  pwagner
! Trying to prevent double-printing of name
!
! Revision 2.94  2009/10/13 00:09:04  pwagner
! Percentages printed with better format; dumpstats.f9h now a direct prerequisite
!
! Revision 2.93  2009/09/10 20:58:00  pwagner
! 3 ways to summarize diffs: 'b' (table), 'r' (rms), 's' (number of differences)
!
! Revision 2.92  2009/08/18 20:41:17  pwagner
! making INTPLACES public can change appearance of dumped ints
!
! Revision 2.91  2009/08/17 16:55:41  pwagner
! Among options string 'd' means 'direct'
!
! Revision 2.90  2009/06/26 00:15:18  pwagner
! Added dumpDumpOptions
!
! Revision 2.89  2009/06/24 22:36:32  pwagner
! Make use of dump1-3.f9h include files
!
! Revision 2.88  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.87  2009/06/16 17:12:55  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.86  2009/05/14 21:59:04  pwagner
! New optional arg AS_GAPS dor dumping 1d bit logicals
!
! Revision 2.85  2009/05/08 00:40:09  pwagner
! Avoid printing blank name when dumped with empty name string
!
! Revision 2.84  2009/04/01 23:30:49  pwagner
! Improved appearance when printing 'nnn lines of xxx not printed.'
!
! Revision 2.83  2008/11/24 19:34:47  pwagner
! Less wasteful of memory; should not segment dault so often
!
! Revision 2.82  2008/10/24 23:21:49  pwagner
! Limits 3d diffs to prevent segment faults
!
! Revision 2.81  2008/08/27 16:23:41  pwagner
! Added dumpSums
!
! Revision 2.80  2008/07/10 00:13:30  pwagner
! SelfDiff can now print half wave lengths
!
! Revision 2.79  2008/07/09 16:30:21  pwagner
! selfdiff can take optional arg unique
!
! Revision 2.78  2008/06/18 20:56:18  pwagner
! New optional arg 'unique' dumps print unique elements, counts only
!
! Revision 2.77  2008/01/09 20:53:22  pwagner
! When NaNs short-circuit dump or diff, print how many and where
!
! Revision 2.76  2008/01/07 21:37:57  pwagner
! Replace DEFAULTUNDEFINEDVALUE with user-settable undefinedValue
!
! Revision 2.75  2007/10/18 23:37:41  pwagner
! Added dumpTable
!
! Revision 2.74  2007/10/12 23:34:49  pwagner
! Using new howfar procedure to summarize diffs after peeling away outliers
!
! Revision 2.73  2007/09/20 18:39:37  pwagner
! Dont interrupt tables of dumped numbers with stamping
!
! Revision 2.72  2007/09/13 21:09:57  pwagner
! -rms and -s combined add new info about & how near
!
! Revision 2.71  2007/07/17 00:21:58  pwagner
! diff stats when rms concentrate on showing ratios
!
! Revision 2.70  2007/07/11 22:29:23  vsnyder
! Add more dumps for complex, add 2x2xN dumps
!
! Revision 2.69  2007/06/14 18:38:36  pwagner
! Allow separate rmsFormat to be set at class level
!
! Revision 2.68  2007/04/14 00:32:47  vsnyder
! Remove OPTIONAL attribute from NAME1 and NAME args of DIFF_...
!
! Revision 2.67  2007/04/03 02:49:23  vsnyder
! Don't look at non-present dummy argument
!
! Revision 2.66  2007/03/23 00:14:30  pwagner
! new optional maxlon arg to dumping n-d chars; sdFormatDefault for numeric dumps public and changeable
!
! Revision 2.65  2007/03/07 21:01:45  pwagner
! Some small changes unrelated to real bugs elsewhere
!
! Revision 2.64  2007/01/31 00:05:43  pwagner
! Added interface for dumping bit arrays
!
! Revision 2.63  2006/08/23 20:06:25  pwagner
! Restored more backward-compatible arglist to DUMP_STRLIST
!
! Revision 2.62  2006/08/22 20:40:04  pwagner
! Added required arg=separator to DUMP_STRLIST
!
! Revision 2.61  2006/07/11 00:24:04  pwagner
! use fillValue properly when computing rms etc.
!
! Revision 2.60  2006/06/24 23:07:04  pwagner
! Changes to reduce memory footprint computing statistics
!
! Revision 2.59  2006/06/09 18:50:12  pwagner
! Avoid dumping an entire array if all elements the same
!
! Revision 2.58  2006/05/24 20:38:14  pwagner
! Allow any of 3 ordering relations for dumping pct
!
! Revision 2.57  2006/04/20 01:09:30  vsnyder
! Don't print lines of zeroes in complex dumps
!
! Revision 2.56  2006/03/22 23:48:52  vsnyder
! Print 'lines ... not printed' instead of 'rows...' to avoid confusion with matrices
!
! Revision 2.55  2006/03/15 17:34:28  pwagner
! Fixed bug causing incorrect rms when diffing with fill values
!
! Revision 2.54  2006/03/03 23:04:55  pwagner
! May dump logical-valued hashes
!
! Revision 2.53  2006/02/28 21:42:29  pwagner
! Replace fillvalues with 0 before computing rms (should actually remove)
!
! Revision 2.52  2006/01/27 01:01:37  pwagner
! Can now dump hashes
!
! Revision 2.51  2005/12/17 00:58:56  pwagner
! dumps of rms, etc. should appear uniform
!
! Revision 2.50  2005/12/16 23:25:13  pwagner
! dumpSize moved from dump0 to output_m
!
! Revision 2.49  2005/12/16 00:04:06  pwagner
! Changes to reflect new MLSFillValues module
!
! Revision 2.48  2005/11/04 18:49:02  pwagner
! Added SelfDiff procedures to dump diffs among consecutive 1d array elems
!
! Revision 2.47  2005/10/03 18:05:52  pwagner
! Allocated memory now dumped in units of MEMORY_UNITS
!
! Revision 2.46  2005/07/20 01:33:47  vsnyder
! Simplify DumpSize routines
!
! Revision 2.45  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.44  2005/05/12 20:43:54  pwagner
! Added diff routines
!
! Revision 2.43  2004/12/14 21:31:59  pwagner
! Optional trim arg added to multi-dim char dumps
!
! Revision 2.42  2004/08/16 17:09:14  pwagner
! 3d integer dumps interface redone like others
!
! Revision 2.41  2004/08/04 23:19:01  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.40  2004/07/23 19:47:20  vsnyder
! Add LBOUND to dump_1d_[double,integer,real]
!
! Revision 2.39  2004/07/23 18:34:59  vsnyder
! Add LBOUND argument to Dump_1d_logical
!
! Revision 2.38  2004/07/21 19:57:21  pwagner
! New rms, wholeArray options to some dumps
!
! Revision 2.37  2004/06/11 19:04:29  pwagner
! which_ints_are_it renamed findAll, moved MLSSets.f90
!
! Revision 2.36  2004/05/27 23:25:33  pwagner
! Added stats parameter; also more 1-d get fillvalue
!
! Revision 2.35  2004/04/05 17:47:43  livesey
! Split dumpsize into real and integer versions
!
! Revision 2.34  2004/04/03 05:43:23  livesey
! Added DumpSize
!
! Revision 2.33  2004/03/30 00:44:10  vsnyder
! Remove unused variable declaration
!
! Revision 2.32  2004/02/26 21:53:31  pwagner
! Can dump ,-separated string list
!
! Revision 2.31  2004/01/21 22:02:19  vsnyder
! Don't use number of entities per line as default width for counter
!
! Revision 2.30  2003/09/19 02:00:14  vsnyder
! More about the goofy Intel compiler
!
! Revision 2.29  2003/09/15 17:43:41  livesey
! Cosmetic change for fussy (and wrong) intel compiler
!
! Revision 2.28  2003/09/06 00:48:40  vsnyder
! Specify default formats with a module parameter instead of literals.
! Change default (1x,1pg13.6) to (1pg14.6) to avoid problems with length
! calculation in output_m.
!
! Revision 2.27  2003/08/08 20:45:42  vsnyder
! Made say_fill_* generic, made them test for numZeroRows, and made them
! optionally do say_subs_only.  This simplified several dump routines.
! Added optional FORMAT arguments in several more routines.
!
! Revision 2.26  2003/07/04 02:41:33  vsnyder
! Substantial simplification by putting little things into subroutines
!
! Revision 2.25  2003/07/02 01:07:27  vsnyder
! Add complex output
!
! Revision 2.24  2003/05/21 19:20:40  vsnyder
! Start a new line after \ 1 if size==1 and clean
!
! Revision 2.23  2003/05/06 00:15:03  pwagner
! Fixed incompatibility with FilterShapes
!
! Revision 2.20.2.4  2003/05/05 23:00:05  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.20.2.3  2003/04/18 20:26:05  vsnyder
! Add WIDTH and FORMAT arguments to 1D_REAL and 1D_DOUBLE
!
! Revision 2.20.2.2  2003/03/27 23:18:33  vsnyder
! Put new-lines in better places
!
! Revision 2.20.2.1  2003/03/14 00:25:47  vsnyder
! Add Dump_2D_Logical, cosmetic changes
!
! Revision 2.20  2002/12/02 23:34:14  pwagner
! Now can dump name/value pairs
!
! Revision 2.19  2002/10/08 00:09:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.18  2002/09/13 18:08:12  pwagner
! May change matrix precision rm from r8
!
! Revision 2.17  2002/02/14 23:21:18  vsnyder
! Work on dumping masks
!
! Revision 2.16  2001/12/08 00:47:51  pwagner
! Added dump_1d_real for s.p. arrays
!
! Revision 2.15  2001/11/29 23:50:53  pwagner
! Added optional blase arg to dump_nd_char; fixed bug where optional
! format not passed from dump_3d_int
!
! Revision 2.14  2001/11/28 23:32:01  livesey
! Fixed bug where dump_2d_integer didn't pass format to 1d dump.
!
! Revision 2.13  2001/10/25 23:30:39  pwagner
! Improved dump_nd_double to skip rows (e.g., of zeros)
!
! Revision 2.12  2001/10/24 18:11:14  pwagner
! which_ints_are_it now works properly
!
! Revision 2.11  2001/10/23 22:40:37  pwagner
! Now dumps 1d,2d,3d char arrays and 3d ints
!
! Revision 2.10  2001/09/28 22:43:20  vsnyder
! Don't print rows of zeroes
!
! Revision 2.9  2001/09/11 22:52:32  livesey
! Added printing of sizes
!
! Revision 2.8  2001/05/11 22:44:54  vsnyder
! Print transpose of 2d-double if it would take fewer lines.  Get rid of
! double printing of "without mask"
!
! Revision 2.7  2001/05/08 20:27:24  vsnyder
! Added an optional 'format' argument in a few more places
!
! Revision 2.6  2001/05/08 17:21:02  livesey
! Added a `clean' option to the array dumps.  This omits the indices at
! the start, making it easier for other programs to read output.
!
! Revision 2.5  2001/05/03 02:12:34  vsnyder
! Insert copyright notice, clean up CVS stuff, cosmetics
!
! Revision 2.4  2001/03/10 03:39:58  vsnyder
! Improve handling of "name" if size==1 or size==0
!
! Revision 2.3  2001/03/02 01:32:08  livesey
! Handles larger arrays better
!
! Revision 2.2  2001/02/28 21:35:27  livesey
! Added dump logical 1d
!
! Revision 2.1  2000/09/13 20:38:50  vsnyder
! Initial code
!
