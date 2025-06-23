! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Dump_1

  ! 1.  Dump routines that use dump_0 but aren't elementary.
  ! 2.  Other dumps that don't use dump_0, and aren't elementary.

  use Dump_Options, only: AfterSub, Clean, Dopts, Laconic, &
    & TheDumpBegins, TheDumpEnds
  use Dump_0, only: Dump, Empty, FinishLine, ILog10, Name_And_Size
  use HighOutput, only: AlignToFit, BlanksToTab, NumToChars, OutputList, &
    & ResetTabs, SetTabs
  use MLSStringLists, only: GetStringElement, NumStringElements
  use MLSStrings, only: LowerCase, WriteIntsToChars
  use Output_m, only: Blanks, NewLine, Output

  implicit none
  private

  public :: Dump, DumpDates, DumpLists, DumpNamedValues, &
    & DumpTable, DumpTextfile

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!     (data types and parameters)

!     (subroutines and functions)
! Dump               (1) Dump a ','-separated string list, one item per line; or
!                    (2) Dump a hash composed of two string lists: keys and values
!                    (3) Dump a hash composed of one string list for keys
!                        and a logical-value array for values
! DumpDates          Dump a 1-d array of tai93 (s. after 1 Jan 1993)
! DumpLists          Dump a 2d array as a set of lists
! DumpNamedValues    Dump an array of paired names and values
! DumpSums           Dump after summing successive array values
!                      ("inverse" of selfDiff)
! DumpTable          Dump a 2-d table of values with headers
! DumpTextfile       Dump the contents of a text file
! === (end of toc) ===

! === (start of api) ===
! DumpDates ( dble dates(:), [int width], [char* dateFormat], &
!     & [char* timeFormat], [log Leapsec] )
! DumpLists ( Array, [char* Name], [int Width], [char* Sep], [char* Delims] )
!       where array can be a 2d array of chars or ints
! dumpNamedValues ( values, strlist names,
!      [char* format, [int width], [char* options] ) 
!       where values can be a 1d array of ints or reals, and
!       names is a string list of corresponding names
! DumpSums ( array, char* name,
!      [fillvalue], [int width], [char* format],
!      [int lbound], [char* options] )
! DumpTable ( values, headers, char* headside
!      [char* format, [char* formats(:)] ) 
!       where values can be a 2d array of ints or reals, and
!       headers is an array the same size as the 2nd index of values
!       format optionally overrides the default format for the numeric type
!       formats allows you to specify a format separately column-by-column
! DumpTextfile ( char* fileName, [char* options] )
!          ..Dump.. interfaces
! Dump_Strlist ( strlist string, char* name, [char* fillvalue],
!      [char* options], [char* InSeparator] )
! Dump_Hash_Log ( log countEmpty, strlist keys, log values[:], char* name, 
!       [char* separator], [char* options] )
! Dump_Hash_Str ( log countEmpty, strlist keys, strlist values, char* name, 
!       [char* separator], [char* options], [char* ExtraValues] )
! === (end of api) ===

  interface Dump
    module procedure Dump_Hash_Log, Dump_Hash_Str, Dump_StrList
  end interface

  interface DumpDates
    module procedure Dump_TAI
  end interface

  interface DumpLists
    module procedure DumpLists_Chars, DumpLists_Ints
  end interface

  interface DumpNamedValues   ! Dump name-value pairs, names in string list
    module procedure DumpNamedValues_Double, DumpNamedValues_Integer
    module procedure DumpNamedValues_Real
  end interface

  interface DumpSums       ! dump after summing successive array values
    module procedure DumpSums_Integer
    module procedure DumpSums_Real
    module procedure DumpSums_Double
  end interface

  interface DumpTable   ! dump table of values, headers
    module procedure DumpTable_Double, DumpTable_Integer
    module procedure DumpTable_Real
  end interface

  ! These are private variables declared module-wide purely for convenience
  integer, parameter :: MaxLineLen              = 120
  integer, parameter :: MaxTextFileLinelen      = 255
  integer, parameter :: MaxTextFileLines        = 600

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! -----------------------------------------------  Dump_Hash_Log  -----
  subroutine Dump_Hash_Log ( CountEmpty, Keys, Values, Name, Separator, Options )
    ! Dumps a hash composed of two string lists: keys and values
    logical, intent(in)          :: CountEmpty
    character(len=*), intent(in) :: Keys
    logical, dimension(:), intent(in) :: Values
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: Separator
    character(len=*), intent(in), optional :: Options

    character( len=len(keys)) :: element
    integer :: J
    integer :: NumElements

    call theDumpBegins ( options )

    NumElements = NumStringElements(keys, dopts(clean)%v, &
     & separator)
    if ( NumElements == 0 ) then
      call empty ( name )
    else
      call name_and_size ( name, .false., NumElements )
      if ( present(name) .and. .not. dopts(laconic)%v ) call newLine
      call output ( '(key)    =   (value)', advance='yes' )
      do j = 1, NumElements
        call GetStringElement( keys, element, j, countEmpty, separator)
        call output ( trim(element), advance='no' )
        call output( ' = ', advance='no' )
        call output ( values(j), advance='no' )
        call finishLine
      end do ! j
    end if
    call theDumpEnds
  end subroutine Dump_Hash_Log

  ! -----------------------------------------------  Dump_Hash_Str  -----
  subroutine Dump_Hash_Str ( CountEmpty, Keys, Values, &
    & Name, Separator, Options, ExtraValues )
    ! Dumps a hash composed of two string lists: keys and values
    logical, intent(in)          :: CountEmpty
    character(len=*), intent(in) :: Keys
    character(len=*), intent(in) :: Values
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: Separator
    character(len=*), intent(in), optional :: Options
    character(len=*), intent(in), optional :: ExtraValues

    character( len=max(len(values), len(keys)) ) :: element
    integer :: J
    integer :: NumElements

    call theDumpBegins ( Options )

    NumElements = NumStringElements(keys, countEmpty, &
     & separator)
    if ( NumElements == 0 ) then
      call empty ( name )
    else
      call name_and_size ( name, .false., NumElements )
      if ( present(name) .and. .not. dopts(laconic)%v ) call newLine
      call output ( '(key)    =   (value)', advance='yes' )
      do j = 1, NumElements
        call GetStringElement( keys, element, j, countEmpty, separator)
        call output ( trim(element), advance='no' )
        call output( ' = ', advance='no' )
        call GetStringElement( values, element, j, countEmpty, separator)
        call output ( trim(element), advance='no' )
        if ( present(ExtraValues) ) then
          call output ( ' : ', advance='no' )
          call GetStringElement( ExtraValues, element, j, countEmpty, separator)
          call output ( trim(element), advance='no' )
        endif
        call finishLine
      end do ! j
    end if
    call theDumpEnds
  end subroutine Dump_Hash_Str

  ! -----------------------------------------------  Dump_StrList  -----
  subroutine Dump_StrList ( String, Name, FillValue, Options, InSeparator )
    ! Dumps a ','-separated string list, one item per line
    ! (Unless it consists of multiple lines)
    character(len=*), intent(in) :: String
    character(len=*), intent(in), optional :: Name
    character(len=*), intent(in), optional :: FillValue
    character(len=*), intent(in), optional :: Options
    character(len=*), optional, intent(in) :: InSeparator

    integer :: J
    integer :: NumElements
    character(len=len(string)) :: myFillValue
    character(len=1), parameter :: CR = ACHAR(13) ! Carriage return
    character(len=1), parameter :: LF = ACHAR(10) ! Line feed
    character(len=1) :: Separator
    logical, parameter :: CountEmpty = .true.
    ! Executable
    if( index(string, CR) > 0 ) then
      call dump( (/ trim(string) /), name, fillvalue, options=options )
      return
    else if( index(string, LF) > 0 ) then
      call dump( (/ trim(string) /), name, fillvalue, options=options )
      return
    end if

    call theDumpBegins ( options )
    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue

    Separator = ','
    if ( present(INSeparator) ) Separator = INSeparator
    
    NumElements = NumStringElements(string, countEmpty, &
     & separator)
    if ( NumElements == 0 ) then
      call empty ( name )
    else if ( NumElements == 1 ) then
      call name_and_size ( name, dopts(clean)%v, 1 )
      call output ( trim(string), advance='no' )
      call finishLine
    else
      call name_and_size ( name, dopts(clean)%v, NumElements )
      if ( present(name) .and. .not. dopts(laconic)%v ) call newLine
      do j = 1, NumElements
        call GetStringElement(string, myFillValue, j, countEmpty, separator)
        call output ( trim(myFillValue), advance='no' )
        call finishLine
      end do ! j
    end if
    call theDumpEnds
  end subroutine Dump_StrList

  ! --------------------------------------------------  DumpLists  -----
  ! Dump a 2d array as a set of lists
  ! Show name (if any)
  subroutine DumpLists_Chars ( Array, Name, Width, Sep, Delims )
    ! Args
    character(len=*), dimension(:,:), intent(in)   :: Array
    character(len=*), intent(in), optional         :: Name
    integer, optional, intent(in)                  :: Width
    character(len=*), optional, intent(in)         :: Sep
    character(len=*), optional, intent(in)         :: Delims
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
    end if
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
    do line = 1, nLines
      n1 = n2 + 1
      n2 = min( n2 + nElemsPerLine, size(array,2) )
      call newLine
      call output ( n1, places=max(4,ilog10(size(array,2)+1)) , advance='no' )
      call output ( afterSub , advance='no' )
      do n=n1, n2
        call blanksToTab
        call outputList( array(:,n), sep, delims)
      end do
    end do
  end subroutine DumpLists_Chars

  subroutine DumpLists_Ints ( ARRAY, Name, WIDTH, SEP, DELIMS )
    ! Args
    integer, dimension(:,:), intent(in)            :: array
    character(len=*), intent(in), optional         :: Name
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
    end if
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
      end do
    end do
    ! Reset tabs
    call newLine
    call resetTabs
  end subroutine DumpLists_Ints

  ! --------------------------------------------  DumpNamedValues  -----
  ! Another hash-like dump:
  ! Show names and related (numerical) values
  ! -------------------------------------  DumpNamedValues_Double  -----
  subroutine DumpNamedValues_Double ( Values, Names, Format, Width, Options )
    double precision, intent(in)           :: Values(:)
    character(len=*), intent(in), optional :: Names
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Width ! How many pairs per line (1)?
    character(len=*), intent(in), optional :: Options
    
    integer :: J, K, L
    integer :: MyWidth
    character(len=24) :: myName
    call theDumpBegins ( options )
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    if ( dopts(clean)%v ) call output(' ', advance='yes')
    l = 0
    do j = 1, size(values), MyWidth
      do k = 1, MyWidth
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

  end subroutine DumpNamedValues_Double

  ! ------------------------------------  DumpNamedValues_Integer  -----
  subroutine DumpNamedValues_Integer ( Values, Names, Format, Width, Options )
    integer, intent(in)                    :: Values(:)
    character(len=*), intent(in), optional :: Names
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Width ! How many pairs per line (1)?
    character(len=*), intent(in), optional :: Options
    
    integer :: J, K, L
    integer :: MyWidth
    character(len=24) :: myName
    call theDumpBegins ( options )
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    if ( dopts(clean)%v ) call output(' ', advance='yes')
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

  end subroutine DumpNamedValues_Integer

  ! ---------------------------------------  DumpNamedValues_Real  -----
  subroutine DumpNamedValues_Real ( Values, Names, Format, Width, Options )
    real, intent(in)                       :: Values(:)
    character(len=*), intent(in), optional :: Names
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Width ! How many pairs per line (1)?
    character(len=*), intent(in), optional :: Options
    
    integer :: J, K, L
    integer :: MyWidth
    character(len=24) :: myName
    call theDumpBegins ( options )
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    if ( dopts(clean)%v ) call output(' ', advance='yes')
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

  end subroutine DumpNamedValues_Real

 ! ----------------------------------------------------  DumpSums  -----
 ! This family of routines dumps the running sum:
 ! summed(i) == ( array(i) + summed(i-1) )
  subroutine DumpSums_Double ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options )
    double precision, intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    double precision, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Lbound ! Low bound for Array
    character(len=*), intent(in), optional :: Options
    ! Local variables
    integer                                  :: i
    double precision, dimension(size(array)) :: summed
    ! Executable
    if ( size(array) < 1 ) return
    summed(1) = array(1)
    if ( size(array) > 1 ) then
      do i=2, size(array)
        summed(i) = array(i) + summed(i-1)
      end do
    end if
    call dump ( summed, Name, &
      & fillValue, width, format, lbound, Options )
  end subroutine DumpSums_Double

  subroutine DumpSums_Integer ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options )
    integer, intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Lbound ! Low bound for Array
    character(len=*), intent(in), optional :: Options
    ! Local variables
    integer                                  :: i
    integer, dimension(size(array)) :: summed
    ! Executable
    if ( size(array) < 1 ) return
    summed(1) = array(1)
    if ( size(array) > 1 ) then
      do i=2, size(array)
        summed(i) = array(i) + summed(i-1)
      end do
    end if
    call dump ( summed, Name, &
      & fillValue, format, width, lbound, Options )
  end subroutine DumpSums_Integer

  subroutine DumpSums_Real ( Array, Name, &
    & FillValue, Width, Format, Lbound, Options )
    real, intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    real, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: Lbound ! Low bound for Array
    character(len=*), intent(in), optional :: Options
    ! Local variables
    integer                                  :: i
    real, dimension(size(array)) :: summed
    ! Executable
    if ( size(array) < 1 ) return
    summed(1) = array(1)
    if ( size(array) > 1 ) then
      do i=2, size(array)
        summed(i) = array(i) + summed(i-1)
      end do
    end if
    call dump ( summed, Name, &
      & fillValue, width, format, lbound, Options )
  end subroutine DumpSums_Real

  ! --------------------------------------------------  DumpTable  -----
  ! This family of routines dumps a table of values and headers
  ! The headside determines on which of the 4 sides 
  ! ('top', 'left', 'right', 'bottom') we place the headers
  subroutine DumpTable_Double ( Values, Headers, HeadSide, Format, Formats )
    double precision, dimension(:,:), intent(in) :: Values
    character(len=*), dimension(:), intent(in)   :: Headers
    character(len=*), intent(in)                 :: HeadSide
    character(len=*), intent(in), optional       :: Format
    character(len=*), dimension(:), intent(in), optional :: Formats
    include 'dumpTable.f9h'
  end subroutine DumpTable_Double

  subroutine DumpTable_Integer ( Values, Headers, HeadSide, Format, Formats )
    integer, dimension(:,:), intent(in)          :: Values
    character(len=*), dimension(:), intent(in)   :: Headers
    character(len=*), intent(in)                 :: HeadSide
    character(len=*), intent(in), optional       :: Format
    character(len=*), dimension(:), intent(in), optional :: Formats
    include 'dumpTable.f9h'
  end subroutine DumpTable_Integer

  subroutine DumpTable_Real ( Values, Headers, HeadSide, Format, Formats )
    real, dimension(:,:), intent(in)             :: Values
    character(len=*), dimension(:), intent(in)   :: Headers
    character(len=*), intent(in)                 :: HeadSide
    character(len=*), intent(in), optional       :: Format
    character(len=*), dimension(:), intent(in), optional :: Formats
    include 'dumpTable.f9h'
  end subroutine DumpTable_Real

  ! ----------------------------------------------------- Dump_TAI -----
  ! Dumps date (re)formatted according to Dates_Module.
  ! This is here instead of in Dates_Module to avoid a circular dependence.
  subroutine Dump_TAI( taiDates, Name, Width, DateFormat, TimeFormat, &
    & LeapSec )
    ! Dump an array of tai dates in whatever formats the user specifies
    ! Warnings  By default:
    ! 
    ! (1) does not correct for leap seconds
    ! (2) we dump both date and time fields
    ! If dateFormat is present and blank or 'none', don't print date
    ! If timeFormat is present and blank or 'none', don't print time
    ! Args
    use Dates_Module, only: MaxUTCStrLength, ReformatDate, ReformatTime, &
      & SplitDateTime, TAI93s2UTC
    use HighOutput, only: OutputNamedValue
    use MLSStrings, only: LowerCase
    double precision, dimension(:)       :: TAIDates ! tai93 (s)
    character(len=*), optional           :: Name
    integer, intent(in), optional        :: Width
    character(len=*), optional           :: DateFormat
    character(len=*), optional           :: TimeFormat
    logical, optional                    :: Leapsec ! Correct for leap seconds?
    ! Internal variables
    character(len=16)                    :: Date, Time
    character(len=maxUTCStrLength), dimension(size(taiDates)) &
      &                                  :: Dates
    logical, parameter                   :: Debug = .false.
    integer                              :: Error ! Should we check on this?
    integer                              :: I
    ! Executable
    if ( size(taidates) < 1 ) then
      call dump ( taiDates, name )
      return
    end if
    do i = 1, size(taiDates)
      dates(i) = tai93s2utc( taiDates(i), leapsec )
    end do
    if ( Debug ) call outputNamedValue( 'last date (utc)', dates(size(taiDates)) )
    if ( present(dateFormat) .or. present(timeFormat) ) then
      do i = 1, size(taiDates)
        call splitDateTime( dates(i), error, date, time )
        if ( present(dateFormat) ) then
          if ( len_trim(dateFormat) < 1 .or. lowercase(dateFormat) == 'none' ) then
            date = ' '
          else
            date = ReformatDate( date, toForm=dateFormat )
          end if
        end if
        if ( present(timeFormat) ) then
          if ( len_trim(timeFormat) < 1 .or. lowercase(timeFormat) == 'none' ) then
            time = ' '
          else
            time = ReformatTime( time, Form=timeFormat )
          end if
        end if
        dates(i) = trim(date) // time
      end do
    end if
    call dump( dates, name, width=width )
  end subroutine Dump_TAI

  ! -----------------------------------------------  DumpTextfile  -----
  ! This dumps a table of values and headers
  ! The headside determines on which of the 4 sides 
  ! ('top', 'left', 'right', 'bottom') we place the headers
  subroutine DumpTextfile ( FileName, Options )
    use HighOutput, only: OutputNamedValue
    use IO_Stuff, only: Read_Textfile
    use Output_m, only: Blanks, Output
    character(len=*), intent(in)           :: FileName
    character(len=*), intent(in), optional :: Options
    ! Internal variables
    character(len=maxTextFileLinelen), dimension(maxTextFileLines) :: Lines
    integer :: K, NLines
    ! Executable
    lines = ' '
    call read_textfile( FileName, lines, nLines=nLines )
    call outputNamedValue ( 'text file name', trim(FileName) )
    call outputNamedValue ( 'number of lines', nLines )
    call blanks( 80, fillChar='-', advance='yes' )
    do k = 1, nLines
      call output ( trim(lines(k)), advance='yes' )
    end do
    call blanks( 80, fillChar='-', advance='yes' )
  end subroutine DumpTextfile

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Dump_1

! $Log$
! Revision 2.9  2019/10/31 22:58:00  pwagner
! Moved DumpException out
!
! Revision 2.8  2019/07/22 22:23:19  pwagner
! Dump_TAI now takes optional arg LeapSec
!
! Revision 2.7  2017/12/07 02:41:07  vsnyder
! Give DEEBUG a value, and make it a parameter
!
! Revision 2.6  2017/07/31 23:07:39  pwagner
! TheDumpEnds now got from dump_options
!
! Revision 2.5  2017/07/10 18:44:44  pwagner
! May print a little more in Dump_TAI
!
! Revision 2.4  2017/01/19 23:31:41  pwagner
! Can DumpException
!
! Revision 2.3  2016/09/22 22:52:42  pwagner
! Dump_Hash_Str can now take optional ExtraValues arg
!
! Revision 2.2  2016/07/28 03:19:52  vsnyder
! Move comments here from dump_0
!
! Revision 2.1  2016/07/28 01:41:48  vsnyder
! Initial Commit
!
