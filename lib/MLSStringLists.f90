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
module MLSStringLists               ! Module to treat string lists
!=============================================================================

  use IO_Stuff, only: PrintMessage
  use MLSCommon, only: BareFNLen, MLSMSG_Error
  use MLSFinds, only: FindFirst, FindLast
  use MLSStrings, only: Capitalize, CompressString, IsAlphabet, LowerCase, &
    & NCopies, ReadIntsFromChars, ReadNumsFromChars, Replace, Reverse, &
    & Reverse_Trim, SplitDetails, Squeeze, StrEq, Trim_Safe, WriteIntsToChars
  use Sort_M, only: Sortp
  implicit none
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

!

! === (start of toc) ===
! This module contains some higher-level string handling stuff for mls
! See below for what we mean by a stringList
! Applications include hdf, the Switches string, and command line arguments

!     c o n t e n t s
!     - - - - - - - -

!     (parameters and data)
! STRINGLISTOPTIONS  Default options string
! KEYNOTFOUND        key not found among keyList
! KEYBEYONDHASHSIZE  index of key in keyList > size(hash array)
! LENORSIZETOOSMALL  Either charsize of strs or size(ints) too small

!     (subroutines and functions)
! Array2List         Converts an array of strings to a single string list
! BooleanValue       Evaluate a boolean formula: e.g. 'p and not {q or r)'
! BuildHash          Builds a hash out of a hash constructor like
!                     [ 100mb : 0.5ppmv, 50mb : 0.1ppmv, 31mb : 0.12ppmv ]
! CapitalizeArray    Capitalize each line of an array of lists
! CapitalizeList     Capitalize the 1st char in each element of a list
! CatLists           cats 2 string lists, taking care if either one is blank
! EvaluateFormula    Evaluates a string formula, substituting character values
!                      for occurrences like ${3} with 3rd arg
! ExpandStringRange  Turns '1,2-5,7' into '1,2,3,4,5,7' or ints or logicals
! ExtractSubString   Extracts portion of string sandwiched between sub1 and sub2
! GetHashElement     Returns value from hash list corresponding to key string
! GetMatchedParens   Returns indexes of matched parens '(' ')'
! GetStringElement   Returns n'th element of string list
! GetUniqueInts      Returns array of only unique entries from input array
! GetUniqueList      Returns str list of only unique entries from input list
! GetUniqueStrings   Returns array of only unique entries from input array
! insertHashElement  Insert a scalar or array-valued element into a hash
! Intersection       Return the intersection of two stringlists; may return blank
! IsInList           Is string in list? options may expand criteria
! List2Array         Converts a single string list to an array of strings
! listMatches        Return list of matches for string in a list
! LoopOverFormula    Looping while it evaluates a string formula, substituting
!                      for occurrences like ${arg} with value of arg
! NCharsInFormat     How many characters in a format spec
! NumStringElements  Returns number of elements in string list
! OptionDetail       Returns detail or arg of option in list of options
! ParseOptions       Parse options from commandline
! PutHashElement     Puts value into hash list corresponding to key string
! ReadIntsFromList   Read an array of ints from a string list
! ReadNumsFromList   Read an array of floats from a string list
! RemoveElemFromList Removes occurrence(s) of elem from a string list
! RemoveHashArray    Removes key strings and corresponding values
!                      based on a named array
! RemoveHashElement  Removes key string and corresponding value
! RemoveListFromList Removes occurrence(s) of elems in a string list from another
! RemoveNumFromList  Removes a numbered elem from a string list
! RemoveOption       Removes an option from a list of options
! RemoveSwitchFromList  
!                    Removes a switch from a list of switches
! ReplaceSubString   Replaces occurrence(s) of sub1 with sub2 in a string
! ReverseList        Turns 'abc,def,ghi' -> 'ghi,def,abc'
! ReverseStrings     Turns (/'abc','def','ghi'/) -> (/'ghi','def','abc'/)
! SnipList           Like RemoveElemFromList, but a function
! SortArray          Turns (/'def','ghi','abc'/) -> (/'abc','def','ghi'/)
! SortList           Turns 'def,ghi,abc' -> 'abc,def,ghi'
! StringElement      Returns string element in string list at index number
! StringElementNum   Returns element number of test string in string list
! SwitchDetail       Returns detail level of switch in list of switches
! Unquote            Removes surrounding [quotes]
! Unwrap             Unwrap a multi-line string to a single line
! Wrap               Wrap a string to fit within prescribed width
!                      using separator as newline
! WriteIntsToList    Write an array of ints as a string list
! === (end of toc) ===

! === (start of api) ===
! Array2List (char* inArray(:), char* outList(:), &
!   & [char inseparator], [int ordering], [char leftRight]) 
! log BooleanValue ( char* str, char* lkeys, log lvalues[:] )
! BuildHash (char* Constructor, char* values, &
!   & [char* operator], [char separator], [char* options])
! char* CatLists (char* str1, char* str2)
! CapitalizeArray ( char* inArray(:), char* outArray(:), &
!   & [char inseparator], [char* ignore] )
! char* CapitalizeList ( char* str, [char inseparator], [char* ignore] )
! ExpandStringRange (char* str, char* outst)
! ExtractSubString (char* str, char* outstr, char* sub1, char* sub2, &
!       & [char* how], [log no_trim])
! char* EvaluateFormula (char* formula, char* values(:), [char* values(:)])
! GetHashElement (hash {keys = values}, char* key, 
!   char* outElement, log countEmpty, [char inseparator], [log part_match])
! GetMatchedParens ( char* str, int pairs(2,:) )
! GetStringElement (strlist inList, char* outElement,
!   nElement, log countEmpty, [char inseparator], [int SeparatorLocation] )
! GetUniqueInts (int ints(:), int outs(:), int noUnique, 
!    [int extra(:)], [int fillValue], [int minValue]) 
! GetUniqueList (char* str, char* outstr(:), int noUnique, &
!   & [char inseparator], [log IgnoreLeadingSpaces], [char* fillValue], &
!   & [char* options])
! GetUniqueStrings (char* inList(:), char* outList(:), int noUnique, 
!   [char* extra(:)], [char* fillValue], [char* options])
! char* intersection ( char* str1, char* str2, [char* options] )
! log IsInList ( strlist stringList, char* string, [char* options] )
! List2Array (strlist inList, char* outArray(:), log countEmpty,
!   [char inseparator], [log IgnoreLeadingSpaces])
! char* listMatches (strlist stringList, char* string, [char* options] )
! LoopOverFormula (char* formula, char* arg, char* values(:), char* results(:))
! int nCharsinFormat ( char* Format )
! int NumStringElements(strlist inList, log countEmpty, &
!   & [char inseparator], [int Longestlen])
! char* optionDetail ( strlist inList, &
!     [char single_option], [char* multi_option], [int pattern], &
!     [char delims(2)] )
! ParseOptions(strlist cmdline, char* opts_out(:), &
!     char single_options(:), int pattern, char* multi_options(:), &
!     [char delims(2)], [char* cmdargs(:)] )
! PutHashElement (hash {keys = values}, char* key, 
!   char* elem, log countEmpty, [char inseparator], [log part_match])
! ReadIntsFromList(strlist inList, int ints(:), &
!    & [char inseparator])
! ReadNumsFromList(strlist inList, num nums(:), &
!    & [char inseparator], [char* ignore], [int error])
! RemoveElemFromList(strlist inList, strlist outList, char* elem, &
!    & [char inseparator], [char* options])
! RemoveHashArray (hash {keys = values}, char* key, 
!   log countEmpty, [char inseparator], [log part_match])
! RemoveHashElement (hash {keys = values}, char* key, 
!   log countEmpty, [char inseparator], [log part_match])
! RemoveListFromList(strlist inList, strlist outList, strlist exclude, &
!    & [char inseparator], [char* options])
! RemoveNumFromList(strlist inList, strlist outList, int nElement, &
!    & [char inseparator], [char* options])
! char* RemoveOption ( strlist inList, strlist outList, char* option, &
!    & [int pattern], [char delims(2)] )
! RemoveSwitchFromList( strlist inList, strlist outList, char* switch, &
!    & [char inseparator], [char* options] )
! ReplaceSubString (char* str, char* outstr, char* sub1, char* sub2, &
!       & [char* which], [log no_trim])
! strlist ReverseList (strlist str, [char inseparator])
! ReverseStrings (char* str[*], char* reverse[*])
! strlist SnipList (strlist str, char* elem, [char inseparator])
! SortArray (char* inStrArray(:), int outIntArray(:), &
!   & [char* sortedArray(:)], [char* options])
! SortList (strlist inStrArray, int outIntArray(:), &
!     [char inseparator], [strlist inStrArray],  [char* options])
! char* StringElement (strlist inList, &
!   nElement, log countEmpty, [char inseparator])
! int StringElementNum(strlist inList, char* test_string, log countEmpty, &
!   & [char inseparator], [log part_match])
! int SwitchDetail(strlist inList, char* test_switch, [char* options])
!      by default, options is "-f" to ignore leading spaces
! char* unquote (char* str, [char* quotes], [char* cquotes], [char* options])
! char* unwrap ( char* str )
! wrap ( char* str, char* outstr, int width, [char inseparator], &
!   & [char break], [char mode], [char* quotes], [int addedLines], &
!   & [log dontSqueeze] )
! WriteIntsToList( int ints(:), strlist List )

! in the above, a string list is a string of elements (usu. comma-separated)
! e.g., units='cm,m,in,ft'
! an array is a Fortran array of strings or integers
! a hash is a list of key strings and their associated values
! (a string list, or array of ints or logicals)
! Many of these routines take optional arguments that greatly modify
! their default operation

! One area of possible improvement, or change, anyway, is the choice of
! commas for separators between elements of a string. This is in accord
! with hdfeos dimension fields, etc. It is not ideal for the most general
! case where, for example, a string element might itself contain a comma.
! In the most general case we ought to allow for, and consider moving
! the default to, a non-ascii character to use for separator, e.g., achar(0)
! or NULL.

! One standard is the character flag "options" which affects how loosely
! string matches may be interpreted, quotes treated, and how string elements are
! counted in lists
! it may include any of the following (poss. in combination, e.g. "-wc")
! w    Wildcard * which allows 'a*' to equal 'abcd'
! b    backward search when that makes sense
! c    case insensitive which allows 'ABCD' to equal 'abcd'
! f    flush left which allows 'abcd' to equal '  abcd'
! e    count consecutive separators as enclosing an empty element
! n    reverse sense of match (where appropriate)
! s{.} use character between braces (here a ".") instead of "," as separator
! k    strict; i.e. remove quotes only if they match
! p    stripany; remove any quotes
! r    remove any quoted sub-strings
! x    extract any quoted substrings
! S    match strings as if Switches; e.g. 'pro' matches 'pro1'
! L    keep Last match in GetUnique instead of first match

! We hope eventually that options will replace the countEmpty, caseSensitive, 
! etc. separate optional args to many of the current module procedures

! An argument can be made that countEmpty should be TRUE by default 
! rather than FALSE
! One example is which FALSE seems best is the slash in file paths
! where we want /data/a.dat and /data//a.dat to be synonyms

! The two procedures EvaluateFormula and LoopOverFormula deserve some
! careful explanation
! We intend them to help pry mlsl2's l2cf away from its chemical dependence
! on m4. Two hurdles to the intervention are m4's statement functions
! and the idea of a loop. Our replacements: EvaluateFormula and LoopOverFormula.
! A typical formula in m4 might look something like
!     T h e   m 4   w a y
!   !define(writeStandardProductWithColumn,{
!   Label, quantity=state.$1, label='$1'
!   Label, quantity=state.column_$1, label='$1 column'
!   DirectWrite, type=l2gp, hdfVersion=!hdfVersion, $
!     file='!l2gpFilename($1)', $
!     source=state.$1, precision=outputPrecision.$1, $
!     status=otherDiagnostics.status$1, quality=otherDiagnostics.quality$1, $
!     source=state.column_$1, precision=outputPrecision.column_$1!nl})
! The equivalent using MLSStringLists would be
!     T h e   M L S S t r i n g L i s t s   w a y
!   character(len=80), dimension(:), parameter :: (/                              &
! "Label, quantity=state.${1}, label='${1}'                                    ", &
! "source=state.${1}, precision=outputPrecision.${1}, $                        ", &
! "status=otherDiagnostics.status${1}, quality=otherDiagnostics.quality${1}, $ ", &
! "source=state.column_${1}, precision=outputPrecision.column_${1}!nl})        "  &
! /)
! You will note that instead of "$1" the argument to be substituted must be
! represented using the "${1}" idiom. In case we were preparing to loop over
! the formula with multiple qtys of the state vector, we would use
! the "${qty}" idiom instead of "${1}".

! Warnings: 
! (1) In the routines Array2List, and SortArray
!     the input arguments include an array of strings;
!     This array is of assumed-size
!     I.e., all elements from array(1:size(array)) are relevant
!     Therefore in calling one of these you probably need to call it as
!       call SortArray(myArray(1:mySize), ..
!     to avoid operating on undefined array elements
! (2) In operating on string lists it is sometimes assumed that no
!     element is longer than a limit: MAXSTRELEMENTLENGTH
! (3) Integer hashes should not be used if some negative
!     values are expected. The value KEYNOTFOUND=-1 is used to indicate
!     "no such key."
! (4) "No such key" is indicated by FALSE for logical values and "," strings
! (5) If the optional extra array or list is supplied to the GetUnique...
!     function, repeated elements purely in the first arg are left undeleted;
!     if you want uniqueness among them, too, you must invoke it twice:
!     first w/o the extra arg, and the second time with the extra arg
! === (end of api) ===

  public :: Array2List, BooleanValue, BuildHash, &
    & CapitalizeArray, CapitalizeList, CatLists, &
    & EvaluateFormula, ExpandStringRange, ExtractSubstring, &
    & GetHashElement, GetMatchedParens, GetStringElement, &
    & GetUniqueInts, GetUniqueStrings, GetUniqueList, &
    & InsertHashElement, Intersection, IsInList, &
    & List2Array, LoopOverFormula, ListMatches, &
    & NCharsInFormat, NumStringElements, &
    & OptionDetail, ParseOptions, PutHashElement, &
    & ReadIntsFromList, ReadNumsFromList, &
    & RemoveElemFromList, RemoveListFromList, RemoveNumFromList, &
    & RemoveHashArray, RemoveHashElement, RemoveOption, RemoveSwitchFromList, &
    & ReplaceSubstring, ReverseList, ReverseStrings, &
    & SnipList, SortArray, SortList, StringElement, StringElementNum, &
    & SwitchDetail, &
    & Unquote, Unwrap, Wrap, WriteIntsToList

! A private type
  type :: Index_Stack_t
    integer :: Index = -1              ! into the expresion
    double precision :: Memory = 0.0d0 ! As accounted in Allocate_Deallocate
    integer :: Sys_Memory = 0          ! In use, in kB (1024), as accounted by
                                       ! the system and accessed by Memory_Used
    integer :: String = 0              ! Index in string table
    integer :: Text = 0                ! Index in string table
    integer :: Tree = 0                ! Where in l2cf, -1 if stack not allocated,
                                       ! -2 if stack index < 1, -3 if stack index
                                       ! > stack_ptr.
  end type
  type(Index_Stack_t), allocatable, save :: Stack(:)
  type(Index_Stack_t)                    :: frame
  integer, save                    :: stack_ptr

  interface BooleanValue
    module procedure BooleanValue_log, BooleanValue_str
  end interface

  interface CatLists
    module procedure CatLists_str, CatLists_int, CatLists_intarray
  end interface

  interface EvaluateFormula
    module procedure EvaluateFormula_string, EvaluateFormula_array
  end interface

  interface GetHashElement
    module procedure GetHashElement_str
    module procedure GetHashElement_strarray
    module procedure GetHashElement_int
    module procedure GetHashElement_log
  end interface

  interface GetStringHashElement
    module procedure GetHashElement_str
    module procedure GetHashElement_strarray
  end interface

  interface ExpandStringRange
    module procedure ExpandStringRange_str, ExpandStringRange_ints, &
      & ExpandStringRange_log, ExpandStringRange_real
  end interface

  interface PutHashElement
    module procedure PutHashElement_str
    module procedure PutHashElement_int
    module procedure PutHashElement_log
    module procedure PutHashElement_strarray
  end interface

  interface MakeStringHashElement
    module procedure PutHashElement_str
    module procedure PutHashElement_strarray
  end interface

  interface ReadNumsFromList
    module procedure ReadRealArrayFromString
    module procedure ReadDoubleArrayFromString
  end interface

  interface RemoveHashElement
    module procedure RemoveHashElement_str
    ! module procedure RemoveHashElement_int
    ! module procedure RemoveHashElement_log
  end interface

  interface Unwrap
    module procedure Unwrap_array, Unwrap_list
  end interface

  interface Wrap
    module procedure Wrap_array, Wrap_sca
  end interface

  ! Public data
  character(len=16), public, save :: STRINGLISTOPTIONS = ' '

  ! Error return values from:
  ! GetHashElement (int args)
  integer, public, parameter      :: KEYNOTFOUND = -1
  integer, public, parameter      :: KEYBEYONDHASHSIZE = KEYNOTFOUND-1
  ! strings2Ints
  integer, public, parameter      :: lenORSIZETOOSMALL = -999
  
  ! A limitation among string list operations
  integer , parameter             :: MAXELEMENTLENGTH    = 80
  integer, private, parameter     :: MaxNumSwitches      = 256
  integer, private, parameter     :: MAXSTRLISTLENGTH    = 4*4096
  integer, private, parameter     :: MAXSTRELEMENTLENGTH = BareFNlen

  character (len=1), parameter    :: COMMA = ','
  character (len=1), parameter    :: BLANK = ' '   ! Returned for any element empty

  logical, private, save          :: countEmpty          
  logical, private, save          :: caseSensitive       
  logical, private, save          :: ignoreLeadingSpaces 
  character(len=1), private, save :: separator           
  logical, parameter              :: deeBug = .false.
  character (len=*), public, parameter :: MLSMSG_Allocate = &
     & "Allocation failed: "
  character (len=*), public, parameter :: MLSMSG_DeAllocate = &
     & "Deallocation failed: "
contains

  ! ---------------------------------------------  Array2List  -----

  ! This subroutine returns a (usually) comma-separated string list, interpreted it
  ! as a list of individual elements, given an equivalent array of
  ! sub-strings in which the n'th element becomes the n'th element

  ! As an optional arg the separator may supplied, in case it isn't a comma
  ! As an optional arg the ordering in which the array elements are to be
  ! taken may be supplied; e.g. (/4, 1, 3, 2/) means 1st take 4th element,
  ! then 1st, then 3rd, and finally 2nd: list[k] = array[ordering[k]]
  ! (unless the further optional arg leftRight is also supplied and equals
  ! one of {"l", "L"} in which case list[ordering[k]] = array[k])

  subroutine Array2List ( inArray, outList, inseparator, ordering, leftRight )
    ! Dummy arguments
    character (len=*), intent(out)                :: outList
    character (len=*), dimension(:), intent(in)   :: inArray
    character (len=*), optional, intent(in)       :: inseparator
    integer, dimension(:), optional, intent(in)   :: ordering
    character (len=1), optional, intent(in)       :: leftRight

    ! Local variables
    integer :: listElem, arrayElem, nElems

    character (len=1)               :: separator
    character (len=1)               :: myLeftRight
    ! Executable code

    if(present(inseparator)) then
      separator = inseparator
    else
      separator = COMMA
    endif

    if(present(leftRight)) then
      myleftRight = Capitalize(leftRight)
    else
      myleftRight = "R"
    endif

    if ( len(outList) <= 0 ) return
    outList = BLANK
    nElems = size(inArray)
    if ( nElems <= 0 ) return
      listElem = 1
    do
      if (.not. present(ordering) ) then
        arrayElem = ListElem
      elseif (myLeftRight == "R") then
        arrayElem = ordering(ListElem)
      else
        ! Try to invert ordering function
        do arrayElem=1, nElems
          if ( ordering(arrayElem) == listElem ) exit
        enddo
        arrayElem = min(arrayElem, nElems)
      endif
      if ( listElem == 1 ) then
        outList = trim(inArray(arrayElem))
      else
        outList = trim(outList) // separator // trim(inArray(arrayElem))
      endif
      listElem = listElem + 1
      if ( listElem > min(nElems, len(outList)) ) return
    enddo

  end subroutine Array2List

  ! ----------------------------------------  BooleanValue  -----
  ! Takes a well-formed formula and returns its value
  ! given a hash of variables and their values
  ! E.g., given 'p and not (q or r)' 
  ! and p=q=TRUE, r=FALSE, returns FALSE
  ! The str and lkeys are all of type character
  ! lvalues may be either an array of logicals or else
  ! a stringlist like lkeys
  ! if a stringlist, it will be converted into an array of logicals
  ! in a natural way
  ! '[tT]*'       -> .true.
  ! anything else -> .false.
  function BooleanValue_log ( str, lkeys, lvalues, separator ) result(BooleanValue)

    ! Method:
    ! Progressively collapse all the '(..)' pairs into their values
    ! until only primitives remain
    ! then evaluate the primitives
    !
    ! Limitations:
    ! does not check for unmatched parens or other illegal syntax
    !
    ! What about precedence? Does it assign higher precedence
    ! to "not" than to "and"? Higher precedence to "and" than to "or"?

    ! Originally it did not. Now it does.

    ! To see what we're talking about
    ! consider the following logical expression
    !    not a and b or c
    ! It should be evaluated the same as
    !    ((not a) and b) or c   <-- right
    ! instead of, e.g.
    !    not (a and (b or c))   <-- wrong
    
    !--------Argument--------!
    character (len=*), intent(in)           :: str
    character (len=*), intent(in)           :: lkeys
    logical, dimension(:), intent(in)       :: lvalues
    logical                                 :: BooleanValue
    character (len=1), optional, intent(in) :: separator
    ! Internal variables
    logical, parameter              :: deeBug = .false.
    integer                         :: level
    integer                         :: level2
    integer, parameter :: MAXNESTINGS=64 ! Max number of '(..)' pairs
    character(len=MAXSTRLISTLENGTH) :: collapsedstr
    integer, dimension(2,1)         :: pairs
    character(len=MAXSTRLISTLENGTH) :: part1
    character(len=MAXSTRLISTLENGTH) :: part2
    character(len=MAXSTRLISTLENGTH) :: part3
    character(len=MAXSTRLISTLENGTH) :: part21
    character(len=MAXSTRLISTLENGTH) :: part22
    character(len=MAXSTRLISTLENGTH) :: part23
    character(len=MAXSTRLISTLENGTH) :: mstr
    logical :: pvalue
    character(len=16) :: vChar

    ! Executable
    BooleanValue = .FALSE.
    if ( str == ' ' ) return
    if ( min(len_trim(lkeys), size(lvalues)) < 1 ) return

    mstr = lowerCase(str)
    ! We're unable to ensure operator precedence
    ! so we'll attempt to identify "and"s and "not"s
    ! and surround such subexpressions with extra parentheses
    
    ! This wasn't sufficient alone,
    ! so we duplicated it inside the endless loop
    ! over nested matched pairs of parentheses.

    call ReorderPrecedence( mstr, collapsedstr )
    if ( DeeBUG ) then
      print *, 'incoming ', trim(mstr)
      print *, 'after reordering precedence ', trim(collapsedstr)
    endif

    ! Collapse every sub-formula nested within parentheses
    do level =1, MAXNESTINGS ! To prevent endlessly looping if ill-formed
      if ( DEEBug ) print *, 'collapsedstr: ', trim(collapsedstr)
      if ( index( collapsedstr, '(' ) < 1 ) exit
      ! call SplitNest ( collapsedstr, part1, part2, part3 )
      call GetMatchedParens( collapsedstr, pairs )
      part1 = ' '
      part2 = ' '
      part3 = ' '
      if ( pairs(1, 1) < 1 ) then
        part2 = collapsedstr
      elseif ( pairs(1, 1) == 1 ) then
        part2 = collapsedstr(2:pairs(2,1)-1)
        if ( pairs(2, 1) < len_trim(collapsedstr) ) &
          & part3 = collapsedstr(pairs(2,1)+1:)
      elseif ( pairs(2, 1) == len_trim(collapsedstr) ) then
        part1 = collapsedstr(1:pairs(1,1)-1)
        part2 = collapsedstr(pairs(1,1)+1:pairs(2,1)-1)
      else
        part1 = collapsedstr(1:pairs(1,1)-1)
        part2 = collapsedstr(pairs(1,1)+1:pairs(2,1)-1)
        part3 = collapsedstr(pairs(2,1)+1:)
      endif
      ! Now evaluate the part2
      if ( DeeBUG ) then
        print *, 'part1 ', trim(part1)
        print *, 'part2 ', trim(part2)
        print *, 'part3 ', trim(part3)
      endif
      ! Hackery-quackery alert:

      ! We're unable to ensure operator precedence
      ! (being too inept or lacking in ideas)
      ! so we'll attempt to identify "and"s and "not"s
      ! and surround such subexpressions with extra parentheses

      call ReorderPrecedence( part2, mstr )
      if ( DEEBug ) then
        print *, 'incoming ', trim(part2)
        print *, 'after reordering precedence ', trim(mstr)
      endif
      
      ! The next block of code is in an endless loop
      ! of its own until there are no more parens in mstr
      ! (which were created in ReorderPrecedence as explained above)
      precedence: do level2 =1, MAXNESTINGS ! To prevent endlessly looping
        if ( DEEBug ) print *, 'mstr: ', trim(mstr)
        if ( index( mstr, '(' ) < 1 ) exit precedence
        ! if ( mstr /= part2 ) then
        call GetMatchedParens( mstr, pairs )
        part21 = ' '
        part22 = ' '
        part23 = ' '
        if ( pairs(1, 1) < 1 ) then
          part22 = mstr
        elseif ( pairs(1, 1) == 1 ) then
          part22 = mstr(2:pairs(2,1)-1)
          if ( pairs(2, 1) < len_trim(mstr) ) &
            & part23 = mstr(pairs(2,1)+1:)
        elseif ( pairs(2, 1) == len_trim(mstr) ) then
          part21 = mstr(1:pairs(1,1)-1)
          part22 = mstr(pairs(1,1)+1:pairs(2,1)-1)
        else
          part21 = mstr(1:pairs(1,1)-1)
          part22 = mstr(pairs(1,1)+1:pairs(2,1)-1)
          part23 = mstr(pairs(2,1)+1:)
        endif
        if ( deeBug ) print *, 'Evaluate part22 primitive ', trim(part22)
        pvalue = evaluatePrimitive( trim(part22) )
        vChar = 'true'
        if ( .not. pvalue ) vChar = 'false'
        if ( DEEBug ) then
          print *, 'part21 ', trim(part21)
          print *, 'part22 ', trim(part22)
          print *, 'part23 ', trim(part23)
          print *, '1st vchar ', trim(vchar)
        endif
        ! And substitute its value for the spaces it occupied
        if (  part21 // part23 == ' ' ) then
          mstr = vChar
        elseif (  part21 == ' ' ) then
          mstr = trim(vChar) // ' ' // part3
        elseif ( part23 == ' ' ) then
          mstr = trim(part21) // ' ' //  vChar
        else
          mstr = trim(part21) // ' ' //  vChar // ' ' // part23
          if ( DEEBug ) then
            print *, 'part21 ', trim(part21)
            print *, '2nd vchar ', trim(vchar)
            print *, 'part23 ', trim(part23)
          endif
        endif
        if ( DEEBug ) then
          print *, 'mstr ', trim(mstr)
        endif
      enddo precedence ! precedence loop
      if ( DEEBug ) then
        print *, 'mstr (after precedence loop)', trim(mstr)
      endif
      part2 = mstr
      if ( part2 == ' ' ) then
        ! This should never happen with well-formed formulas
        collapsedstr = part1
        cycle
      else
        pvalue = evaluatePrimitive( trim(part2) )
        vChar = 'true'
        if ( .not. pvalue ) vChar = 'false'
      endif
      ! And substitute its value for the spaces it occupied
      if (  part1 == ' ' ) then
        collapsedstr = trim(vChar) // ' ' // part3
      elseif ( part3 == ' ' ) then
        collapsedstr = trim(part1) // ' ' // vChar
      else
        collapsedstr = trim(part1) // ' ' // trim(vChar) // &
          & ' ' // part3
      endif
      if ( DeeBUG ) then
        print *, 'collapsedstr ', trim(collapsedstr)
      endif
    enddo
    ! Presumably we have collapsed all the nested '(..)' pairs by now
    BooleanValue = evaluatePrimitive( trim(collapsedstr) )
  contains
  function isBalanced ( str ) result ( itIs )
    ! Are numbers of '(' and ')' in str equal?
    ! See also Enclosure subroutine in MLSStrings
    character(len=*), intent(in) :: str
    logical                      :: itIs
    ! Internal variables
    integer :: i
    integer :: balance
    ! Executable
    ! 1st--Just count
    ! Not as strict as "do they balance", but easier to check
    itIs = ncopies( trim(str), '('           ) ==  &
      &    ncopies( trim(str), ')'           ) 
    if ( .not. itIs ) return
    ! Now look more closely
    ! Each '(' adds 1, each ')' subtracts 1, though never going negative
    ! If 0 when done, we're balanced, != 0 if not
    balance = 0
    do i=1, len_trim(str)
      select case(str(i:i))
      case ('(')
        balance = balance + 1
      case (')')
        balance = max( 0, balance - 1 )
      ! case default (leaves balance unchanged)
      end select
    enddo
    itIs = ( balance == 0 )
  end function isBalanced

  function isParenthetical ( str ) result ( itIs )
    ! TRUE if 1st non-blank is '(' and last non-blank is ')'
    ! Do we really need to make this a recursive function?
    character(len=*), intent(in) :: str
    logical                      :: itIs
    itIs = index( adjustl(str), '('           ) == 1 .and. &
      &    index( trim(str), ')', back=.true. ) == len_trim(str)
  end function isParenthetical

    subroutine ReorderPrecedence ( arg, sult )
      character(len=*), intent(in)                  :: arg
      character(len=*), intent(out)                 :: sult
      !
      integer                                       :: i
      integer                                       :: n
      character                                     :: separator
      character(len=MAXSTRLISTLENGTH)               :: element
      character(len=MAXSTRLISTLENGTH)               :: rev
      character(len=MAXSTRLISTLENGTH)               :: tmp
      !
      sult = ' '
      if ( len_trim(arg) < 1 ) then
        return
      endif
      separator = '+'
      ! 1st, replace each instance of ' or ' with '+'
      !      replace each instance of ' and ' with '*'
      call ReplaceSubString( arg, rev, ' or ', ' + ', &
        & which='all', no_trim=.true. )
      call ReplaceSubString( rev, tmp, ' and ', ' * ', &
        & which='all', no_trim=.true. )
      if ( DeeBUG ) print *, 'arg in ReorderPrecedence: ', trim(arg)
      n = NumStringElements( tmp, COUNTEMPTY, inseparator='+' )
      do i=1, n
        call GetStringElement ( tmp, element, i, countEmpty, inseparator='+' )
        ! Be careful not to surround a string that ends with an op
        ! (Why does this keep happening?)
        rev = Reverse_trim(element)
        ! Is it an op?
        if ( index( '+-*/', rev(1:1) ) > 0 ) then
          sult = catLists( sult, element, inseparator='+' )
          cycle
        endif
        ! Surround term with parentheses if it's a product or quotient or '^'
        ! but not if it's (already) parenthetical or parentheses are not balanced
        if ( ( index(element, '*') > 0  ) .and. &
          & ( &
          &   .not. isParenthetical(element) .and. isBalanced(element) &
          & ) &
          & ) then
          if ( DEEBug ) then
          print *, 'element to be parenthesized: ', trim(element)
          print *, 'isParenthetical(element): ', isParenthetical(element)
          print *, 'isBalanced(element): ', isBalanced(element)
          print *, 'ncopies( trim(element), "("): ', ncopies( trim(element), '(')
          print *, 'ncopies( trim(element), ")"): ', ncopies( trim(element), ')')
          endif
          element = '(' // trim(element) // ')'
        endif
        ! Avoid producing two consecutive ops, like '+ +' or '* +'
        ! which would happen if collapsedstr ends with an op
        ! Begin by finding its last non-blank char
        !
        rev = Reverse_trim(element)
        ! Is it an op?
        if ( index( '+-*/', rev(1:1) ) > 0 ) then
          sult = trim(sult) // element
        else
          sult = catLists( sult, element, inseparator='+' )
        endif
      enddo
      ! Now replace each '+' with ' or '
      ! and each '*' with ' and '
      call ReplaceSubString( sult, rev, '+', ' or ', &
        & which='all', no_trim=.true. )
      call ReplaceSubString( rev, sult, '*', ' and ', &
        & which='all', no_trim=.true. )
    end subroutine ReorderPrecedence

    function evaluatePrimitive(primitive) result(value)
      ! Evaluate an expression composed entirely of
      ! (0) constants (e.g., 'true')
      ! (1) primitives (e.g., 'p')
      ! (2) unary operators ('not')
      ! (3) binary operators ('or', 'and', 'xor')
      ! Dummy args
      character(len=*) :: primitive
      logical          :: value
      ! Internal variables
      logical          :: done
      integer          :: elem
      logical          :: hit
      character(len=3) :: lastOp ! "or", "and", "xor"
      integer          :: n
      logical          :: negating
      logical          :: part
      character(len=32) :: variable
      ! Executable
      value = .true.
      done = .false.
      negating = .false.
      elem = 0
      lastOp = 'nul' ! 'or'
      n = NumStringElements( trim(primitive), countEmpty=.false., inseparator=' ' )
      do
        ! go through the elements, re-evaluating every time we "hit" a primitive
        ! Otherwise revising our lastOp or negating status
        elem = elem + 1
        call GetStringElement ( trim(primitive), variable, elem, &
          & countEmpty=.false., inseparator=' ' )
        select case(trim(variable))
        case ('t', 'true')
          part = .true.
          hit = .true.
        case ('f', 'false')
          part = .false.
          hit = .true.
        case ('or')
          lastOp = 'or'
          hit = .false.
        case ('xor')
          lastOp = 'xor'
          hit = .false.
        case ('and')
          lastOp = 'and'
          hit = .false.
        case ('not', '~')
          negating = .true.
          hit = .false.
        case default
          call GetHashElement( lkeys, lvalues, trim(variable), part, &
            & countEmpty=.true., inseparator=separator )
          hit = .true.
          if ( DeeBug ) then
            print *, trim(variable), ' is ', part
            if ( present(separator) ) print *, 'iachar(separator)', ' is ', iachar(separator)
            ! print *, 'keys', ' is ', trim(lkeys)
            ! print *, 'values', ' is ', lvalues
          endif
        end select
        if ( hit ) then
          if ( negating ) part = .not. part
          select case(lastOp)
          case ('nul')
              value = part
          case ('or')
              value = value .or. part
          case ('and')
              value = value .and. part
          case ('xor')
              value = ( value .or. part ) .and. &
                & .not. ( value .and. part )
          case default
            ! How could this happen?
              call PrintMessage( MLSMSG_Error, ModuleName, &
                & lastOp // ' not a legal binary op in evaluatePrimitive' )
          end select
          negating = .false.
        endif
        if ( DeeBUG ) then
          print *, 'variable ', variable
          print *, 'part ', part
          if ( hit ) print *, 'hit ', hit
          if ( negating ) print *, 'negating ', negating
          print *, 'value ', value
        endif
        done = ( elem >= n )
        if ( done ) exit
      enddo
    end function evaluatePrimitive
  end function BooleanValue_log

  function BooleanValue_str ( str, lkeys, strvalues, separator ) result(BooleanValue)
    !--------Argument--------!
    character (len=*), intent(in)           :: str
    character (len=*), intent(in)           :: lkeys
    character (len=*), intent(in)           :: strvalues
    logical                                 :: BooleanValue
    character (len=1), optional, intent(in) :: separator
    ! Internal variables
    logical, parameter                :: countEmpty = .true.
    logical, dimension(:), pointer    :: lvalues
    integer :: key
    integer :: nkeys
    integer :: status
    ! Executable
    BooleanValue = .false.
    if ( len_trim(str) < 1 ) return
    nkeys = NumStringElements( lkeys, countEmpty, inseparator=separator )
    if ( nkeys < 1 ) return
    nullify( lvalues )
    allocate( lvalues(nkeys), stat=status )
    if ( status /= 0 ) call PrintMessage( MLSMSG_Error, ModuleName, &
      & 'Unable to allocate lvalues in BooleanValue_str' )
    do key = 1, nkeys
      ! print *, 'value(key) ', &
      !   & adjustl(lowercase(StringElement ( strvalues, key, countEmpty, inseparator=separator )))
      lvalues(key) = index( &
        & adjustl(lowercase(&
        & StringElement ( strvalues, key, countEmpty, inseparator=separator )&
        & )) &
        & , 't') == 1
    enddo
    BooleanValue = BooleanValue_log ( str, lkeys, lvalues, separator )
    deallocate ( lvalues, stat=status )
    if ( status /= 0 ) call PrintMessage( MLSMSG_Error, ModuleName, &
      & 'Unable to deallocate lvalues in BooleanValue_str' )
  end function BooleanValue_str

  ! -------------------------------------------------  BuildHash  -----
  ! Build a hash from a constructor
  ! E.g., given '[ 100mb : 0.5ppmv, 50mb : 0.1ppmv, 31mb : 0.12ppmv ]'
  ! Returns 
  ! keys = '100mb, 50mb, 31mb'  
  ! hash = '0.5ppmv, 0.1ppmv, 0.12ppmv'
  ! operator   optionally use something other than ':'
  ! separator  optionally use something other than ','
  ! options    'a' append to existing keys, values
  ! Use:
  ! The constructor example shown above is used in level 2 profile Fills
  subroutine BuildHash( Constructor, Keys, values, operator, separator, options )
    ! Dummy arguments
    character (len=*), intent(in)             :: Constructor
    character (len=*), intent(inout)          :: Keys
    character (len=*), intent(inout)          :: values
    character (len=*), optional, intent(in)   :: operator ! defaults to ':'
    character (len=*), optional, intent(in)   :: separator ! defaults to ','
    character (len=*), optional, intent(in)   :: options
    ! Internal variables
    logical :: append
    logical, parameter :: countEmpty = .false.
    integer :: i
    character(len=128) :: istr
    character(len=3) :: op
    character(len=1) :: sep
    character(len=MAXSTRLISTLENGTH) :: str
    ! Executable
    op = ':'
    if ( present(operator) ) op = operator
    sep = ','
    if ( present(separator) ) sep = separator
    append = .false.
    if ( present(options) ) append = ( index(options, 'a') > 0 )
    if ( .not. append ) then
      keys = ' '
      values = ' '
    endif
    str = unquote( Constructor, quotes='[]$', options='-p' )
    do i=1, NumStringElements( str, countEmpty, inseparator=trim(sep) )
      call GetStringElement( str, istr, i, countEmpty, inseparator=trim(sep) )
      if ( len_trim(istr) < 1 ) cycle
      keys = CatLists( keys, stringElement( istr, 1, countEmpty, inseparator=trim(op) ) )
      values = CatLists( values, stringElement( istr, 2, countEmpty, inseparator=trim(op) ) )
    enddo
  end subroutine BuildHash

  ! -------------------------------------------------  CapitalizeArray  -----
  subroutine CapitalizeArray ( inArray, outArray, inseparator, ignore )
    ! Capitalize the 1st char of line of an array
    ! If size(inseparator) > 1, do the same for each of its elements
    !--------Argument--------!
    character (len=*), dimension(:), intent(in)            :: inArray
    character (len=*), dimension(:), intent(out)           :: outArray
    character (len=*), optional, intent(in)                :: ignore ! Don't change any of these
    character (len=1), dimension(:), optional, intent(in)  :: inseparator
    ! Local variables
    integer                                                :: i
    integer                                                :: j
    ! Executable
    if ( .not. present(inseparator) ) then
      do i=1, size(inArray)
        outArray(i) = CapitalizeList( inArray(i), ignore=ignore )
      enddo
      return
    endif
    do i=1, size(inArray)
      outArray(i) = CapitalizeList( inArray(i), inseparator(1), ignore )
    enddo
    if ( size(inseparator) < 2 ) return
    do j=2, size(inseparator)
      do i=1, size(inArray)
        outArray(i) = CapitalizeList( outArray(i), inseparator(j), ignore )
      enddo
    enddo
  end subroutine CapitalizeArray

  ! -------------------------------------------------  CapitalizeList  -----
  function CapitalizeList ( STR, inseparator, ignore ) result ( OUTSTR )
    ! Capitalize the 1st char of each element of a list
    ! E.g., given str = 'hot,dog,2tommy' 
    ! returns 'Hot,Dog,2Tommy'
    !--------Argument--------!
    character (len=*), intent(in)                 :: STR
    character (len=1), optional, intent(in)       :: inseparator
    character (len=*), optional, intent(in)       :: ignore ! Don't change any of these
    character (len=len(str)+9) :: OUTSTR

    !----------Local vars----------!
    integer                             :: i
    integer                             :: j
    integer                             :: m
    integer                             :: n
    character (len=len(str)+9)          :: temp
    character (len=1)                   :: separator
    character, dimension(len(str)+9)    :: temp1d
    !----------executable part----------!
    if(present(inseparator)) then
      separator = inseparator
    else
      separator = comma
    endif
    n = NumStringElements ( str, countEmpty, inseparator )
    outstr = ' '
    if ( n < 1 ) return
    do i=1, n
      call GetStringElement ( str, temp, i, countEmpty, inSeparator )
      m = len_trim(temp)
      temp1d(1:m) = temp(1:m)
      j = FindFirst ( IsAlphabet(temp1d(1:m)) )
      if ( IsInList ( ignore, temp, options='-fc' ) ) then
        ! Nothing to do here
      elseif ( j > 0 ) then
        temp(j:j) = Capitalize(temp(j:j))
      endif
      outstr = CatLists( outstr, trim_safe(temp), inseparator=inseparator )
    enddo
  end function CapitalizeList

  ! -------------------------------------------------  CatLists  -----
  ! This family lets you Build up a String List from scratch
  ! Typical usage
  ! List = ' ' ! Intialize String List
  ! do i = 1, n
  !   List = CatLists ( List, str(i) ) ! Cat str(i) ono List end
  ! enddo
  !
  ! -------------------------------------------------  CatLists_int  -----
  function CatLists_int (STR1, INT, inseparator) result (OUTSTR)
    ! appends an int onto end of a string list, taking care if it is blank
    ! E.g., given str1 = 'a,b,c' and int = 4
    ! returns 'a,b,c,4'
    ! if str1 is blank, returns just '4'
    !--------Argument--------!
    character (len=*), intent(in) :: STR1
    integer, intent(in)           :: INT
    character (len=*), optional, intent(in)       :: inseparator
    character (len=len(str1)+9) :: OUTSTR

    !----------Local vars----------!
    character (len=1)                :: separator
    character (len=32), dimension(1) :: str2
    !----------executable part----------!
    if(present(inseparator)) then
      separator = inseparator
    else
      separator = comma
    endif
    call writeIntsToChars( (/int/), str2 )
    str2(1) = adjustl(str2(1))
    if ( len_trim(str2(1)) < 1 ) then
      outstr=str1
    elseif ( len_trim(str1) < 1 ) then
      outstr=str2(1)
    else
      outstr = trim(str1) // separator // trim(str2(1))
    endif
  end function CatLists_int

  ! ---------------------------------------------  CatLists_intarray  -----
  function CatLists_intarray (STR1, INTS, inseparator) result (OUTSTR)
    ! appends array of ints onto end of a string list, 
    ! taking care if it is blank
    ! E.g., given str1 = 'a,b,c' and ints = (/4,5,5,0/)
    ! returns 'a,b,c,4,5,5,0'
    ! if str1 is blank, returns just '4,5,5,0'
    !--------Argument--------!
    character (len=*), intent(in)                 :: STR1
    integer, dimension(:), intent(in)             :: INTS
    character (len=*), optional, intent(in)       :: inseparator
    character (len=len(str1)+8*size(ints))        :: OUTSTR, TMPSTR

    !----------Local vars----------!
    integer :: i
    !----------executable part----------!
    outstr = str1
    if ( size(ints) < 1 ) then
      return
    endif
    do i = 1, size(ints)
      tmpstr = outstr
      outstr = CatLists_int(tmpstr, ints(i))
    enddo
  end function CatLists_intarray

  ! -------------------------------------------------  CatLists_str  -----
  function CatLists_str (STR1, STR2, inseparator) result (OUTSTR)
    ! cats 2 string lists, taking care if either is blank
    ! E.g., given str1 = 'a,b,c' and str2 = 'd,e,f'
    ! returns 'a,b,c,d,e,f'
    ! if either is blank, returns the other
    ! if both blank, returns a blank
    !
    !--------Argument--------!
    character (len=*), intent(in) :: STR1
    character (len=*), intent(in) :: STR2
    character (len=*), optional, intent(in)       :: inseparator
    character (len=len(str1)+len(str2)+1) :: OUTSTR

    !----------Local vars----------!
    character (len=1)               :: separator
    !----------executable part----------!
    if(present(inseparator)) then
      separator = inseparator
    else
      separator = comma
    endif
    if ( len_trim(str2) < 1 ) then
      outstr=str1
    elseif ( len_trim(str1) < 1 ) then
      outstr=str2
    else
      outstr = trim(str1) // separator // trim(str2)
    endif
  end function CatLists_str

  ! -------------------------------------------------  EvaluateFormula  -----
  ! Evaluates a string formula, plugging in the nth value for
  ! each occurrence of the nth arg appearing as '${n}'
  ! E.g., if strFun is 
  ! "x${1}: vector, template=state"
  ! and the first arg is "InitPtan" then outStr will be
  ! "xInitPtan: vector, template=state"
  ! 
  ! If keys is present, then instead of '${n}' substitute for
  ! each occurrence of the nth arg appearing as '${key(n)}'
  function EvaluateFormula_array ( FORMULA, VALUES, KEYS ) result (OUTSTR)
    !--------Argument--------!
    character (len=*), intent(in)                         :: FORMULA
    character (len=*), dimension(:), intent(in)           :: VALUES
    character (len=*), dimension(:), optional, intent(in) :: KEYS
    character (len=len(formula)+size(values)*len(values)) :: OUTSTR

    !----------Local vars----------!
    integer :: i, n
    character (len=len(formula)+size(values)*len(values)) :: tmpstr
    character (len=len(values))                           :: value       
    character(len=128)                                    :: variable    
    !----------executable part----------!
    n = size(values)
    outstr = formula
    ! Check whether we have any work to do
    if ( n < 1 .or. index( formula, '${' ) < 1 ) return
    do i=1, n
      value = values(i)
      if ( present(keys) ) then
        variable = '${' // trim(keys(i)) // '}'
      else
        call WriteIntsToChars( i, variable )
        variable = '${' // trim(adjustl(variable)) // '}'
      endif
      tmpstr = outstr
      call ReplaceSubString( tmpStr, outstr, trim(variable), trim(value), &
        & which='all', no_trim=.true. )
    enddo
  end function EvaluateFormula_array

  function EvaluateFormula_string ( formula, values, keys, separator ) &
    & result (outstr)
    !--------Argument--------!
    character (len=*), intent(in)                     :: formula
    character (len=*), intent(in)                     :: values
    character (len=*), optional, intent(in)           :: keys
    character (len=1), optional, intent(in)           :: separator
    character (len=len(formula)+len(values))          :: outstr

    !----------Local vars----------!
    integer :: status, n
    character (len(values)), dimension(:), pointer    :: valuesArray              
    character (len(values)), dimension(:), pointer    :: keysArray                
    logical, parameter                                :: countEmpty = .true.
    !----------executable part----------!
    outstr = formula
    ! Check whether we have any work to do
    if ( len_trim(values) < 1 .or. index( formula, '${' ) < 1 ) return
    n = NumStringElements( values, countEmpty, separator )
    allocate ( valuesArray(n), STAT=status )
    call list2Array( values, valuesArray, countEmpty, inseparator=separator )
    if ( present(keys) ) then
      allocate ( keysArray(n), STAT=status )
      call list2Array( keys, keysArray, countEmpty, inseparator=separator )
      outstr = EvaluateFormula_array ( formula, valuesArray, keysArray )
      deallocate( keysArray, stat=status )
    else
      outstr = EvaluateFormula_array ( formula, valuesArray )
    endif
    deallocate( valuesArray, stat=status )
  end function EvaluateFormula_string

  ! ----------------------------------------  ExpandStringRange_ints  -----
  subroutine ExpandStringRange_ints (instr, ints, LENGTH)
    ! Takes a range and returns an array of integers
    ! E.g., given '1,2-5,7' returns (/1,2,3,4,5,7/)
    !--------Argument--------!
    character (len=*), intent(in) :: instr
    integer, dimension(:), intent(out) :: ints
    integer, optional, intent(out) :: LENGTH  ! number of ints returned
    ! Internal variables
    integer :: elem
    integer :: ErrTyp
    integer :: nelem
    character(len=MAXSTRLISTLENGTH) :: expandedstr
    character(len=16) :: iChar
    logical, parameter :: countEmpty=.true.
    ! Executable
    ints = 0
    if ( present(LENGTH) ) LENGTH = 0
    if ( min(len_trim(instr), size(ints)) < 1 ) return
    call ExpandStringRange_str (instr, expandedstr)
    nelem = NumStringElements(trim(expandedstr), countEmpty)
    if ( nelem < 1 ) return
    do elem = 1, min(nelem, size(ints))
      call GetStringElement (trim(expandedstr), iChar, elem, countEmpty)
      read(iChar, *, iostat=ErrTyp) ints(elem)
    enddo
    if ( present(LENGTH) ) LENGTH = min(nelem, size(ints))
  end subroutine ExpandStringRange_ints

  ! ----------------------------------------  ExpandStringRange_log  -----
  subroutine ExpandStringRange_log (instr, logs, fits, highfit, sense)
    ! Takes a range and returns an array of logicals true for indices
    ! fitting the range (unless sense is false)
    ! E.g., given '1,2-5,7' returns (/T,T,T,T,T,F,T/) (because '6' is outside)
    !--------Argument--------!
    character (len=*), intent(in) :: instr
    logical, dimension(:), intent(out) :: logs
    integer, optional, intent(out) :: fits  ! number of fits in range (trues)
    integer, optional, intent(out) :: highfit  ! index of highest fit
    logical, optional, intent(in) :: sense   ! if false, return F for fits
    ! Internal variables
    integer :: elem
    integer :: ErrTyp
    integer :: indx
    integer :: nelem
    integer :: nfits
    integer :: hfit
    character(len=MAXSTRLISTLENGTH) :: expandedstr
    character(len=16) :: iChar
    logical, parameter :: countEmpty=.true.
    logical :: mySense
    ! Executable
    mySense = .true.
    if ( present(sense) ) mySense = sense
    logs = .not. mySense ! .false.
    if ( present(fits) ) fits = 0
    if ( present(highfit) ) highfit = 0
    if ( min(len_trim(instr), size(logs)) < 1 ) return
    call ExpandStringRange_str (instr, expandedstr)
    nelem = NumStringElements(trim(expandedstr), countEmpty)
    if ( nelem < 1 ) return
    nfits = 0
    hfit = 0
    do elem = 1, min(nelem, size(logs))
      call GetStringElement (trim(expandedstr), iChar, elem, countEmpty)
      read(iChar, *, iostat=ErrTyp) indx
      if ( indx > 0 .and. indx <= size(logs) ) then
        if ( logs(indx) .neqv. mySense ) nfits = nfits + 1
        logs(indx) = mySense ! .true.
        hfit = max(hfit, indx)
      endif
    enddo
    if ( present(fits) ) fits = nfits
    if ( present(highfit) ) highfit = hfit
  end subroutine ExpandStringRange_log

  ! ----------------------------------------  ExpandStringRange_real  -----
  subroutine ExpandStringRange_real (instr, reals, LENGTH)
    ! Takes a range and returns an array of reals
    ! E.g., given '1,2.5-3.5+0.5,7' returns (/1.0,2.5,3.0,3.5,7.0/)
    !--------Argument--------!
    character (len=*), intent(in) :: instr
    real, dimension(:), intent(out) :: reals
    integer, optional, intent(out) :: LENGTH  ! number of reals returned
    ! Internal variables
    integer :: elem
    integer :: ErrTyp
    integer :: nelem
    character(len=MAXSTRLISTLENGTH) :: expandedstr
    character(len=16) :: iChar
    logical, parameter :: countEmpty=.true.
    ! Executable
    reals = -999.99
    if ( present(LENGTH) ) LENGTH = 0
    if ( min(len_trim(instr), size(reals)) < 1 ) return
    call ExpandStringRange_str (instr, expandedstr)
    nelem = NumStringElements(trim(expandedstr), countEmpty)
    if ( nelem < 1 ) return
    do elem = 1, min(nelem, size(reals))
      call GetStringElement (trim(expandedstr), iChar, elem, countEmpty)
      read(iChar, *, iostat=ErrTyp) reals(elem)
    enddo
    if ( present(LENGTH) ) LENGTH = min(nelem, size(reals))
  end subroutine ExpandStringRange_real

  ! --------------------------------------------------  ExpandStringRange_str  -----
  subroutine ExpandStringRange_str (instr, outstr)
    ! Takes a string list and expands any ranges found within:
    ! Ranges are marked by patterns
    ! (1a--simple) 'n-m' integers from n to m inclusive, where m > n
    ! (1b--stride) 'n-m+s' integers from n to m with stride s, where s > 0
    ! (2a--simple) 'n-m' reals from n to m inclusive, where m > n
    ! (2b--stride) 'n-m+s' reals from n to m with stride s, where s > 0
    ! Examples:
    ! '1-10' becomes '1,2,3,4,5,6,7,8,9,10'
    ! '1-5+0.5' becomes '1,1.5,2,2.5,3,3.5,4,4.5,5'
    ! '2-21+2' bcecomes '2,4,6,8,10,12,14,16,18,20' (missing '21')

    ! Errors and limitations:
    ! '5-5' becomes simply '5' (range is inclusive, but not duplicative)
    ! '6-4' becomes '' (we can't go backward)
    ! m, n must be non-negative
    ! The separator between elements must be ','
    ! Although countEmpty set to .true. below, actual behavior
    ! ignores blank elements; e.g. 
    ! '0,1,,2,,4-6,,' becomes '0,1,2,4,5,6'
    ! '1.e-3' type notation confuses the range finder, so we will
    ! go through and replace 
    ! 'e-' -> 'em'
    ! 'e+' -> 'ep'
    ! Args
    character (len=*), intent(in) :: instr
    character (len=*), intent(inout) :: outstr
    ! Internal variables
    character (len=len(outstr)) :: str
    character (len=len(outstr)) :: tempstr
    integer :: dashpos
    integer :: elem
    character(len=*), parameter :: em     = 'em'
    character(len=*), parameter :: eminus = 'e-'
    character(len=*), parameter :: ep     = 'ep'
    character(len=*), parameter :: eplus  = 'e+'
    integer :: ErrTyp
    integer :: m                    ! if substring is of form 'n-m'
    character (len=16) :: mChar
    integer :: n                    ! substring is '..,n[-..],'
    character (len=16) :: nChar
    integer :: nelem
    integer :: pluspos
    integer :: s                    ! if substring is of form 'n-m+s'
    character (len=16) :: sChar
    integer :: t
    character (len=16) :: tChar
    ! These are the real counterparts of n, m s
    real :: rm, rn, rs
    integer :: ns
    logical, parameter :: countEmpty=.true.
    ! Executable
    outstr = instr
    nelem = NumStringElements(instr, countEmpty)
    ! Try to deal with '1.e-(+)..' problem
    str = lowercase(instr)
    call ReplaceSubString( str, tempstr, eminus, em, which='all' )
    call ReplaceSubString( tempstr, str, eminus, em, which='all' )
    dashpos = index(str, '-')
    if ( nelem < 1 .or. dashpos < 1 ) return
    outstr = ' '
    do elem = 1, nelem
      call GetStringElement (str, tempstr, elem, countEmpty)
      dashpos = index(trim(tempstr), '-')
      if ( dashpos < 1 ) then
        ! Must reverse any '1.e-(+)' substitutions we might have made
        call ReplaceSubString( tempstr, nChar, ep, eplus, which='first' )
        call ReplaceSubString( nChar, tempstr, em, eminus, which='first' )
      else
        ! The range operator '-' is invoked
        n = -999
        m = -999
        s = 1
        call GetStringElement (trim(tempstr), nChar, 1, countEmpty, &
          & inseparator='-')
        pluspos = index(trim(tempstr), '+')
        ! print *, 'pluspos: ', pluspos
        if ( pluspos > 0 ) then
          call ExtractSubString (trim(tempstr), mChar, '-', '+')
          ! print *, 'mChar: ', mChar
          call GetStringElement (trim(tempstr), sChar, 2, countEmpty, &
            & inseparator='+')
          ! print *, 'sChar: ', sChar
          ! non-simple pattern is 'n-m+s'
        else
          ! simple pattern is 'n-m'
          call GetStringElement (trim(tempstr), mChar, 2, countEmpty, &
            & inseparator='-')
          sChar = ''
        endif
        ! print *, 'nChar: ', trim(nChar)
        ! print *, 'mChar: ', trim(mChar)
        ! print *, 'sChar: ', trim(sChar)
        if ( index(tempstr, '.') < 1 ) then
          ! print *, 'Case (1a) or (1b)'
          read(nChar, *, iostat=ErrTyp) n
          if ( ErrTyp == 0 ) read(mChar, *, iostat=ErrTyp) m
          if ( ErrTyp == 0 .and. sChar /= '') read(sChar, *, iostat=ErrTyp) s
          ! print *, 'n, m, s: ', n, m, s
          tempstr = ''
          if ( m >= n ) then
            do t=n, m, s
              write(tChar, '(i16)') t
              ! print *, 'tChar: ', trim(tChar)
              tempstr = CatLists(trim(tempstr), adjustl(tChar))
            enddo
          endif
        else
          ! print *, 'Case (2a) or (2b)'
          ! Must reverse any '1.e-(+)' substitutions we might have made
          call ReplaceSubString( nChar, tempstr, em, eminus, which='first' )
          call ReplaceSubString( tempstr, nChar, ep, eplus, which='first' )
          read(nChar, *, iostat=ErrTyp) rn
          call ReplaceSubString( mChar, tempstr, em, eminus, which='first' )
          call ReplaceSubString( tempstr, mChar, ep, eplus, which='first' )
          if ( ErrTyp == 0 ) read(mChar, *, iostat=ErrTyp) rm
          rs = 1.0
          call ReplaceSubString( sChar, tempstr, em, eminus, which='first' )
          call ReplaceSubString( tempstr, sChar, ep, eplus, which='first' )
          if ( ErrTyp == 0 .and. sChar /= '') read(sChar, *, iostat=ErrTyp) rs
          ! print *, 'n, m, s: ', rn, rm, rs
          ns = (rm - rn) / rs + 1.1 ! To deal with roundoff
          tempstr = ''
          do s=1, ns
            write(tChar, '(g10.3)') rn + (s-1)*rs
            ! print *, 'tChar: ', trim(tChar)
            tempstr = CatLists(trim(tempstr), adjustl(tChar))
          enddo
        endif
      endif
      outstr = CatLists(trim(outstr), adjustl(tempstr))
    enddo
  end subroutine ExpandStringRange_str

  ! --------------------------------------------------  ExtractSubString  -----
  subroutine ExtractSubString (instr, outstr, sub1, sub2, how, no_trim)
    ! Takes a string and extracts what is sandwiched between sub1 and sub2
    ! Defaults to choosing only the first occurrence of sub1 and sub2
    ! But if how == 'greedy' chooses last occurrence of sub2
    ! or if how == 'stingy' chooses last occurrence of sub1
    ! Note that, depending on how, we extract:
    !    (let sub1='abc' sub2='def' str='abcabc123defdef')
    ! (a) if how == default => 'abc123'
    ! (b) if how == greedy => 'abc123def'
    ! (c) if how == stingy => '123'
    ! if no_trim is TRUE, sub1 and sub2 may have trailing spaces
    ! that will not be trimmed before attempting to match
    ! Method:
    ! Replace substrings sub1 and sub2 with separator character
    ! and then use GetStringElement to get subelement number 2
    ! We are careful to choose as separator one that is not already present
    ! in the string
    !  
    ! Notes and limitations:
    ! A fundamental issue arises if sub2 occurs before sub1 in the string
    ! do we want to interpret the request such that we
    ! (1) return a blank
    ! (2) look for occurrences of sub2 in the string after sub1
    ! I think we should aim for 2, as it produces a generalization
    ! of picking elements out of a comma-separated list
    !
    ! Misc questions
    ! (1) Will this still work if sub1 has leading or trailing blanks? 
    ! (2) How about sub2?
    ! (3) do we need an optional arg, no_trim, say, that will leave them?
    !     Tried coding it, but can't say for sure it works
    ! (4) What if sub1 is a substring of sub2, or vice versa?
    ! (5) Should we switch to non-ascii characters for use as separator?
    !--------Argument--------!
    character (len=*), intent(in) :: instr
    character (len=*), intent(in) :: sub1
    character (len=*), intent(in) :: sub2
    character (len=*), intent(INOUT) :: outstr
    character (len=*), intent(in), optional :: how
    logical, intent(in), optional :: no_trim

    !----------Local vars----------!
    character (len=len(instr)) :: str
    integer, parameter         :: EARLYSUB2INTERPRETATION = 2 ! 2 or 1
    integer :: i, isub1, isub2, strlen, tmpstrlen
    character (len=7) :: my_how
    character(len=1) :: separator
    character(len=*), parameter :: separators =',.$%#{}()'
    character (len=max(len(instr), len(outstr))) :: tmpstr
    character (len=max(len(instr), len(outstr))) :: tmpstr2
    logical :: my_no_trim, trimming
    !----------Executable part----------!
    my_how = 'first'
    if ( present(how) ) my_how = lowercase(how)
    my_no_trim = .false.
    if ( present(no_trim) ) my_no_trim = no_trim
    trimming = .not. my_no_trim
    outstr = ' '
    strlen = len_trim(instr)
    if (strlen < 1 .or. instr == ' ') return
    ! Which interpretation of sub2 occurring before sub1 do we make?
    isub1 = index(instr, trim(sub1))
    isub2 = index(instr, trim(sub2))
    if ( isub2 < isub1 ) then
      if ( isub2 == 0 ) then
        return
      elseif ( EARLYSUB2INTERPRETATION == 2 ) then
        ! zap every occurrence of sub2 up to position isub1
        call ReplaceSubString (instr(1:isub1-1), tmpstr, sub2, '', &
          & which='all', no_trim=.false.)
        tmpstrlen = len_trim(tmpstr)
        tmpstrlen = len(trim(tmpstr))
        str = ''
        if ( tmpstrlen < 1 ) then
          str = instr(isub1:strlen)
        else
          str = tmpstr(1:tmpstrlen) // instr(isub1:strlen)
        endif
      else
        str = instr
      endif
    else
      str = instr
    endif
    if ( trimming ) then
      if (len_trim(sub1) < 1 &
        & .or. &
        & len_trim(sub2) < 1 .or. index(str, trim(sub1)) == 0 &
        & .or. &
        & index(str, trim(sub2)) == 0 ) RETURN
    else
      if (index(str, sub1) == 0 &
        & .or. &
        & index(str, sub2) == 0 ) RETURN
    endif
    do i=1, len(separators)
      if ( index(str, separators(i:i)) == 0 ) exit
    enddo
    if ( i > len(separators) ) return   ! This means our method will fail
    separator = separators(i:i)
    select case (trim(my_how))
    case ('greedy')
      if ( trimming ) then
        call ReplaceSubString (Reverse(trim(str)), tmpstr, &
          & Reverse(trim(sub2)), separator)
        tmpstr2 = Reverse(trim(tmpstr))
        call ReplaceSubString (tmpstr2, tmpstr, sub1, separator)
      else
        call ReplaceSubString (Reverse(str), tmpstr, &
          & Reverse(sub2), separator, no_trim=.true.)
        tmpstr2 = Reverse(tmpstr)
        call ReplaceSubString (tmpstr2, tmpstr, sub1, separator, no_trim=.true.)
      endif
    case ('stingy')
      if ( trimming ) then
        call ReplaceSubString (Reverse(trim(str)), tmpstr, &
          & Reverse(trim(sub1)), separator)
        tmpstr2 = Reverse(trim(tmpstr))
        call ReplaceSubString (tmpstr2, tmpstr, sub2, separator)
      else
        call ReplaceSubString (Reverse(str), tmpstr, &
          & Reverse(sub1), separator, no_trim=.true.)
        tmpstr2 = Reverse(tmpstr)
        call ReplaceSubString (tmpstr2, tmpstr, sub2, separator, no_trim=.true.)
      endif
    case default
      call ReplaceSubString (str, tmpstr2, sub1, separator, &
        & which='first', no_trim=no_trim)
      call ReplaceSubString (tmpstr2, tmpstr, sub2, separator, &
        & which='first', no_trim=no_trim)
    end select
    call GetStringElement (tmpstr, outstr, 2, .true., &
      & inseparator=separator )

  end subroutine ExtractSubString

  ! ---------------------------------------------  GetStringElement  -----

  ! This subroutine takes a (usually) comma-separated string list,
  ! interprets it as a list of individual elements and returns the
  ! sub-string which is the n'th element
  ! if n is too large or small, it returns the separator
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-separated lists

  ! if countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! if TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the separator may supplied, in case it isn't comma.
  ! Its length can be more than one, in which case any element of it is
  ! a separator.  For example, it could be ", " to find string elements
  ! separated either by commas or spaces.

  ! See also SplitWords

  subroutine GetStringElement ( inList, outElement, nElement, countEmpty, &
    & inseparator, SeparatorLocation )
    ! Dummy arguments
    character (len=*), intent(in)   :: inList
    character (len=*), intent(out)  :: outElement
    integer, intent(in)             :: nElement ! Entry number to return
    logical, intent(in)             :: countEmpty
    character (len=*), optional, intent(in)       :: inseparator
    integer, optional, intent(out)  :: SeparatorLocation ! -1 if no element

    ! Local variables
    integer :: i           ! Loop counters
    integer :: elem, nextseparator

    character (len=1)               :: separator
    ! Executable code

    if(present(inseparator)) then
      separator = inseparator
    else
      separator = COMMA
    endif

    if(nElement.LE.0) then
      outElement = separator
      return
    elseif(len(inList) < nElement) then
      outElement = separator
      return
    endif
    i = 1
    elem = 1
    if ( present(separatorLocation) ) separatorLocation = -1
    do
      if ( i > len(inList) ) then
        outElement = separator
        return
      endif
      nextseparator = i - 1 + SCAN(inList(i:), separator)

      ! No more separators
      if(nextseparator == i - 1) then
        if(elem >= nElement) then
          outElement = inList(i:)
          if ( present(separatorLocation) ) separatorLocation = i
        else
          outElement = separator
        endif
        RETURN

        ! Next separator is the adjacent char
      elseif(nextseparator == i) then
        if(countEmpty) then
          if(elem >= nElement) then
            outElement = BLANK
            RETURN
          else
            elem = elem+1
          endif
        endif

        ! Until next separator is the next element
        else
          if(elem >= nElement) then
            if(i < nextseparator) then
              outElement = inList(i:nextseparator-1)
              if ( present(separatorLocation) ) separatorLocation = i
            else
              outElement = separator
            endif
            RETURN
          elseif(nextseparator >= len(inList)) then
            outElement = separator
            RETURN
          else
            elem = elem+1
          endif
        endif
        i = nextseparator+1
      enddo

  end subroutine GetStringElement

  ! ---------------------------------------------  GetHashElement  -----

  ! This family of subroutines interpret two arguments as a set
  ! of {key = value} pairs

  ! Takes two (usually) comma-separated string lists, interprets
  ! each as a list of elements, treating the first as keys and the second as
  ! a hash table, associative array or dictionary
  ! It returns the sub-string from the hash table corresponding to the key
  ! if the key is not found in the array of keys, it returns the separator
  
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-separated lists

  ! if countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! if TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the separator may supplied, in case it isn't comma
  ! Another optional arg, part_match, returns a match for the 
  ! first hash element merely found in the key; e.g.
  ! 'won, to, tree' and key 'protocol.dat' matches 'to'

  ! Basic premise: Use StringElementNum on key in keyList to find index
  ! Use this index to GetStringElement from HashList

  ! Someday you may wish to define a StringHash_T made up of the two
  ! strings
  
  subroutine GetHashElement_int( KEYS, VALUES, KEY, VALUE, &
  & COUNTEMPTY, INSEPARATOR, PART_MATCH )

  ! if no match found, return KEYNOTFOUND

  ! Someday you may wish to define a StringHash_T made up of the two
  ! strings
  
    ! Dummy arguments
    character (len=*), intent(in)             :: keys
    integer, dimension(:), intent(in)         :: values
    character (len=*), intent(in)             :: key
    integer, intent(out)                      :: value
    logical, intent(in)                       :: countEmpty
    character (len=*), optional, intent(in)       :: inseparator
    logical, optional, intent(in)             :: part_match

    ! Local variables
    integer :: elem

    ! Executable code

    value = KEYNOTFOUND
    elem = StringElementNum(keys, key, countEmpty, inseparator, part_match)
    if( elem <= 0 ) then
    elseif( elem > size(values) ) then
      value = KEYBEYONDHASHSIZE
    else
      value = values(elem)
    endif

  end subroutine GetHashElement_int

  subroutine GetHashElement_log( KEYS, VALUES, KEY, VALUE, &
  & COUNTEMPTY, INSEPARATOR, PART_MATCH )

  ! if no match found, return FALSE

  ! Someday you may wish to define a StringHash_T made up of the two
  ! strings
  
    ! Dummy arguments
    character (len=*), intent(in)             :: keys
    logical, dimension(:), intent(in)         :: values
    character (len=*), intent(in)             :: key
    logical, intent(out)                      :: value
    logical, intent(in)                       :: countEmpty
    character (len=*), optional, intent(in)       :: inseparator
    logical, optional, intent(in)             :: part_match

    ! Local variables
    integer :: elem

    ! Executable code

    value = .FALSE.
    elem = StringElementNum(keys, key, countEmpty, inseparator, part_match)
    if( elem <= 0 ) then
    elseif( elem > size(values) ) then
      value = .false.
    else
      value = values(elem)
    endif

  end subroutine GetHashElement_log

  subroutine GetHashElement_str( KEYLIST, HASHLIST, KEY, OUTELEMENT, &
  & COUNTEMPTY, INSEPARATOR, PART_MATCH )
  
    ! Dummy arguments
    character (len=*), intent(in)   :: keyList
    character (len=*), intent(in)   :: hashList
    character (len=*), intent(in)   :: key
    character (len=*), intent(out)  :: outElement
    logical, intent(in)   :: countEmpty
    character (len=*), optional, intent(in)       :: inseparator
    logical, optional, intent(in)             :: part_match

    ! Local variables
    integer :: elem
    character (len=1)                          :: separator

    ! Executable code

    if(present(inseparator)) then
      separator = inseparator
    else
      separator = COMMA
    endif

    elem = StringElementNum(keyList, key, countEmpty, inseparator, part_match)
    if(elem <= 0) then
      outElement = separator
    else
      CALL GetStringElement(hashList, outElement, elem, &
        & countEmpty, inseparator)
    endif

  end subroutine GetHashElement_str

  subroutine GetHashElement_strarray( KEYLIST, HASHLIST, KEY, ARRAY, &
  & COUNTEMPTY, INSEPARATOR, PART_MATCH )
    ! We fill an array of values from a hash
    ! assuming the array is stored like this
    ! array "name" contains "value_1", "value_2" .. "nalue_n"
    !    key        value
    !    ---        -----
    ! "namen"        "n" (number of elements)
    ! "name(1)"   "value_1"
    ! "name(2)"   "value_2"
    !    .    .    .
    ! "name(n)"   "value_n"
    ! Dummy arguments
    character (len=*), intent(in)                 :: KEYLIST
    character (len=*), intent(in)                 :: HASHLIST
    character (len=*), intent(in)                 :: KEY
    character (len=*), dimension(:), intent(out)  :: ARRAY
    logical, intent(in)                           :: COUNTEMPTY
    character (len=*), optional, intent(in)       :: INSEPARATOR
    logical, optional, intent(in)                 :: PART_MATCH

    ! Local variables
    integer                                       :: j
    character (len=16)                            :: keyString
    integer                                       :: n
    character (len=8)                             :: nCh
    character (len=1)                             :: separator

    ! Executable code

    if(present(inseparator)) then
      separator = inseparator
    else
      separator = COMMA
    endif
    
    array = separator
    keyString = trim(key) // 'n'
    call GetHashElement_str(keyList, hashList, keyString, nCh, &
      & countEmpty, inseparator, part_match)
    if ( nCh == separator ) return
    call readIntsFromChars ( nCh, n )
    
    do j=1, n
      call writeIntsToChars( j, nCh )
      keyString = trim(key) // '(' // trim(adjustl(nCh)) // ')'
      call GetHashElement_str( keyList, hashList, keyString, array(j), &
        & countEmpty, inseparator, part_match )
    enddo

  end subroutine GetHashElement_strarray

! -------------------------------------------------  GetMatchedParens  -----
  subroutine GetMatchedParens ( str, pairs, numpairs )
    ! Get the indexes of Matched Parens in str
    ! sorted so that the most deeply nested comes first.
    ! If shape(pairs) = (/2, 1/)
    ! then return only that most deeply nested pair.
    ! Args:
    character(len=*), intent(in)              :: str
    integer, dimension(:,:), intent(out)      :: pairs
    integer, optional      , intent(out)      :: numpairs
    ! Internal variables
    integer                                   :: i
    integer                                   :: n
    integer                                   :: k ! shape(pairs) = (2,k)
    type(Index_Stack_t)                       :: frame
    ! Executable
    k = size(pairs, 2)
    pairs = 0
    n = 0
    if ( present(numpairs) ) numpairs = 0
    do i=1, len_trim(str)
      if ( str(i:i) == '(' ) then
        call Push ( i )
      elseif ( str(i:i) == ')' ) then
        call Pop ( frame )
        ! call outputnamedvalue ( 'matched parens', (/ frame%index, i /) )
        ! call outputnamedvalue ( 'sub-str', &
        !  & str(frame%index+1:i-1) )
        n = n + 1
        pairs(1, n) = frame%index
        pairs(2, n) = i
        if ( present(numpairs) ) numpairs = n
        if ( n >= k ) exit
      endif
    enddo
    call Deallocate_Index_Stack
  end subroutine GetMatchedParens

  ! ---------------------------------------------  GetUniqueInts  -----

  ! This subroutine takes an array of ints and returns another containing
  ! only the unique entries. The resulting array is supplied by the caller
  ! Its first noUnique entries return the unique values found
  ! Later entries are returned unchanged
  ! if optional extra array is supplied, instead
  ! returns entries from first array not also found in second
  ! if optional fillValue is supplied, values = fillValue are ignored
  ! else if optional minValue is supplied, values < minValue are ignored
  ! Some checking is done to make sure it's appropriate

  subroutine GetUniqueInts(ints, outs, noUnique, extra, fillValue, minValue)
    ! Dummy arguments
    integer, dimension(:) :: ints
    integer, dimension(:) :: outs
    integer :: nounique ! number of unique entries
    integer, optional, dimension(:) :: extra
    integer, optional               :: fillValue
    integer, optional               :: minValue

    ! Local variables
    integer :: i,j,k           ! Loop counters
    logical, dimension(:), allocatable :: duplicate ! Set if already found
    integer :: status        ! Status from allocate

    integer :: extrasize
    integer :: howmanymax
    integer :: insize

    ! Executable code, setup arrays
    inSize=SIZE(ints)
    allocate (duplicate(inSize), STAT=status)
    if (status /= 0) CALL PrintMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"duplicate in GetUniqueInts")
    if ( present(extra) ) then
      extraSize=size(extra)
      howManyMax = inSize
      ! print *, 'SIZE(inList) ', inSize
      ! print *, 'SIZE(extra) ', extraSize
    else
      extraSize = -1
      howManyMax = inSize-1 ! don't bother with last one
    endif
    duplicate = .FALSE.

    ! Go through and find duplicates

    do i = 1, howManyMax
       if (.NOT. duplicate(i)) then
         if ( extraSize < 1 ) then
          do j = i+1, inSize
             if (ints(j)==ints(i)) duplicate(j)=.TRUE.
          end do
         else
          do j = 1, extraSize
             if (extra(j)==ints(i)) duplicate(i)=.TRUE.
          end do
         endif
       endif
    end do

    ! Ignore any values = fillValue
    if ( present(fillValue) ) then
      duplicate = duplicate .or. (ints == fillValue)
    elseif ( present(minValue) ) then
      duplicate = duplicate .or. (ints < minValue)
    endif
    ! Count how many unique ones there are

    noUnique=count(.NOT. duplicate)

    if (noUnique>SIZE(outs)) CALL PrintMessage(MLSMSG_Error,ModuleName, &
         & "outs too small in GetUniqueInts")

    if ( noUnique > 0 ) then
      j=1
      UniqueLoop: do i = 1, noUnique
         k = findFirst(.not. duplicate(j:))
         ! print *, 'j: ', j, '   k: ', k
         if ( k+j-1 > inSize ) then
           call PrintMessage(MLSMSG_Error, ModuleName, &
             & "k goes past array end in GetUniqueInts")
           outs(i)=ints(inSize)
           return
         elseif ( k > 0 ) then
           outs(i)=ints(k+j-1)  ! was ints(j)
           j = j + k
         else
           exit UniqueLoop
         endif
         ! j=j+1
         if ( j > inSize ) exit UniqueLoop
      end do UniqueLoop
    endif

    deallocate ( duplicate )
  end subroutine GetUniqueInts

  ! ---------------------------------------------  GetUniqueList  -----

  ! This subroutine takes a string list and returns another containing
  ! only the unique entries. The resulting list is supplied by the caller
  ! (You may safely use the same variable for str and outStr) (?? f95 std??)
  ! E.g., given 'one,two,three,one,four' returns 'one,two,three,four'
  ! if optional string list str2 is supplied, instead
  ! returns list from str that are not also in str2
  ! if optional FillValue supplied, ignores any entries = fillvalue

  subroutine GetUniqueList( str, outStr, noUnique, &
    & inseparator, IgnoreLeadingSpaces, str2, fillValue, options )
    ! Dummy arguments
    character (len=*), intent(in) :: str
    character (len=*), intent(out) :: outstr
    integer :: nounique ! number of unique entries
    ! logical, intent(in)                           :: countEmpty
    character (len=*), optional, intent(in)       :: inseparator
    logical, optional, intent(in)       :: IgnoreLeadingSpaces
    character (len=*), optional, intent(in)       :: str2
    character (len=*), optional, intent(in)       :: fillValue
    character (len=*), optional, intent(in)       :: options

    ! Local variables
    logical :: countEmpty
    character(len=8) :: myOptions
    character (len=MAXSTRELEMENTLENGTH), dimension(:), allocatable    &
      &                             :: inStringArray, outStringArray, inStrAr2
    integer :: nElems
    integer :: nElems2
    integer :: Longestlen
    integer :: status

    ! Executable code
    myOptions = ' '
    if ( present(options) ) myOptions = options
    countEmpty = ( index(myOptions, 'e') > 0 ) 
    if ( len(str) <= 0 .or. len(outstr) <= 0 ) return
    nElems = NumStringElements(str, countEmpty, inseparator, Longestlen)
    noUnique = nElems
    if ( present(str2) ) then
      outStr = ''
      if ( nElems < 1 ) return
    else
      outStr = str
      if ( nElems <= 1 ) return
    endif
    if ( Longestlen > MAXSTRELEMENTLENGTH ) then
      ! print *, 'str: ', trim(str)
      ! print *, 'len(str): ', len(str)
      ! print *, 'Longestlen: ', Longestlen
      ! print *, 'nElems: ', nElems
      call PrintMessage(MLSMSG_Error, ModuleName, &
         & "Element LENGTH too long in GetUniqueList")
      return
    endif
    allocate (inStringArray(nElems), outStringArray(nElems), STAT=status)
    ! print *, 'shapes: ', &
    !   & (/ size(inStringArray), size(outStringArray) /)
    if (status /= 0) CALL PrintMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"stringArray in GetUniqueList")
    call list2Array(str, inStringArray, countEmpty, inseparator, &
     & IgnoreLeadingSpaces)
    ! print *, 'str ', str
    ! print *, 'nElems ', nElems
    ! do status=1, nElems
    !   print *,  inStringArray( status )
    ! enddo
    if ( present(str2) ) then
      nElems2 = NumStringElements(str2, countEmpty, inseparator, Longestlen)
      allocate (inStrAr2(nElems2), STAT=status)
      if (status /= 0) CALL PrintMessage(MLSMSG_Error,ModuleName, &
           & MLSMSG_Allocate//"stringArray2 in GetUniqueList")
      call list2Array(str2, inStrAr2, countEmpty, inseparator, &
       & IgnoreLeadingSpaces)
      call GetUniqueStrings( inStringArray, outStringArray, noUnique, &
       & inStrAr2, fillValue, options )
      if ( noUnique > 0 ) then
        call Array2List(outStringArray(1:noUnique), outStr, &
         & inseparator)
      else
        outStr=''
      endif
      deallocate (inStringArray, outStringArray, inStrAr2)
    else
      ! print *, 'About to getUniqueStrings'
      call GetUniqueStrings( inStringArray, outStringArray, noUnique, &
      & fillValue=fillValue, options=options )
      ! print *, 'noUnique: ', noUnique
      ! do status=1, noUnique
      !   print *,  outStringArray( status )
      ! enddo
      if ( noUnique > 0 ) then
        call Array2List(outStringArray(1:noUnique), outStr, &
         & inseparator)
      else
        outStr=''
      endif
      deallocate (inStringArray, outStringArray)
    endif
      ! print *, 'done with getUniqueList'
  end subroutine GetUniqueList

  ! ---------------------------------------------  GetUniqueStrings  -----

  ! This subroutine takes an array of strings and returns another containing
  ! only the unique entries. The resulting array is supplied by the caller
  ! if optional extra array is supplied, instead
  ! returns entries from first array not also found in second
  ! if optional FillValue supplied, ignores any entries = fillvalue
  ! Some checking is done to make sure it's appropriate

  subroutine GetUniqueStrings( inList, outList, noUnique, &
    & extra, fillValue, options )
    ! Dummy arguments
    character (len=*), dimension(:) :: inList
    character (len=*), dimension(:) :: outList
    integer :: noUnique ! Number of unique entries
    character (len=*), optional, dimension(:) :: extra
    character (len=*), optional, intent(in)       :: fillValue
    character (len=*), optional, intent(in)       :: options

    ! Local variables
    logical, dimension(:), allocatable :: duplicate ! Set if already found

    integer :: extraSize
    integer :: howManyMax
    integer :: i,j,k           ! Loop counters
    integer :: inSize
    logical :: keepLast
    character(len=len(inList)), dimension(size(inList)) :: list
    character(len=8) :: myOptions
    integer :: status        ! Status from allocate
    logical :: Switchable

    ! Executable code, setup arrays
    ! print *, 'Now ino getUniqueStrings'
    myOptions = ' '
    if ( present(options) ) myOptions = options
    Switchable = ( index(myOptions, 'S') > 0 )
    keepLast = ( index(myOptions, 'L') > 0 )
    inSize=SIZE(inList)
    allocate (duplicate(inSize), STAT=status)
    if (status /= 0) CALL PrintMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"duplicate")
    if ( present(extra) ) then
      extraSize=size(extra)
      howManyMax = inSize
      ! print *, 'SIZE(inList) ', inSize
      ! print *, 'SIZE(extra) ', extraSize
    else
      extraSize = -1
      howManyMax = inSize-1 ! don't bother with last one
    endif
    duplicate = .FALSE.

    ! if we are keeping the last instance of a duplicated string, we'll
    ! simply reverse the order of the list, keep the first instance of that
    ! reversed list, then reverse again at the end
    if ( keepLast ) then
      call reverseStrings( inList, list )
    else
      list = inList
    endif
    ! Go through and find duplicates
    do i = 1, howManyMax
       if (.NOT. duplicate(i)) then
         if ( extraSize < 1 ) then
          do j = i+1, inSize
             ! if (List(j)==List(i)) duplicate(j)=.TRUE.
             duplicate(j) = duplicate(j) .or. matchem( List(j), List(i) )
             ! if ( duplicate(j) ) print *, List(j), List(i)
          end do
         else
          do j = 1, extraSize
             ! if (extra(j)==List(i)) duplicate(i)=.TRUE.
             duplicate(j) = duplicate(j) .or. matchem( extra(j), List(i) )
          end do
         endif
       endif
       ! print *, i, duplicate(i)
    end do

    ! Ignore any values = fillValue
    if ( present(fillValue) ) then
      ! duplicate = duplicate .or. (List == fillValue)
    do i = 1, size(duplicate)
      duplicate(i) = duplicate(i) .or. matchem( List(i), fillValue )
    enddo
    endif
    ! Count how many unique ones there are

    noUnique=count(.NOT. duplicate)

    if (noUnique>SIZE(outList)) CALL PrintMessage(MLSMSG_Error,ModuleName, &
         & "outList too small")
    if (len(outList)<len(List)) CALL PrintMessage(MLSMSG_Error,ModuleName, &
         & "outList strings too small")
    outList=""

    ! print *, 'NoUnique: ', noUnique
    if ( noUnique > 0 ) then
      ! do j=1, inSize, 20
      !   print *, (duplicate(j+i), i=0, min(19, inSize-j))
      ! enddo
      j=1
      UniqueLoop: do i = 1, noUnique
         ! UniqueHuntLoop: do
         !   if (.NOT. duplicate(j)) EXIT UniqueHuntLoop
         !   j=j+1
         !   if ( j > inSize ) exit UniqueLoop
         ! end do UniqueHuntLoop
         k = findFirst(.not. duplicate(j:))
         ! print *, 'j: ', j, '   k: ', k
         if ( k+j-1 > inSize ) then
           call PrintMessage(MLSMSG_Error, ModuleName, &
             & "k goes past array end in GetUniqueStrings")
           outList(i)=List(inSize)
           return
         elseif ( k > 0 ) then
           outList(i)=List(k+j-1)  ! was List(j)
           j = j + k
         else
           exit UniqueLoop
         endif
         ! j=j+1
         if ( j > inSize ) exit UniqueLoop
      end do UniqueLoop
    endif
    ! print *, 'done with UniqueLoop'
    ! if we reversed the order, recover the original order
    if ( keepLast ) then
      list = outList
      ! print *, 'reversing strings'
      call reverseStrings( list(1:noUnique), outList(1:noUnique) )
    endif

    deallocate ( duplicate )
    ! print *, 'Leaving getUniqueStrings'
    contains
    function matchem( str1, str2 ) result ( match )
      ! Test for match between str1 and str2 according to options
      ! Args
      character(len=*), intent(in) :: str1
      character(len=*), intent(in) :: str2
      logical                      :: match
      ! Internal variables
      character(len=len(str1)) :: switch1, switch2
      integer :: details1, details2
      ! Executable
      match = .false.
      if ( .not. Switchable ) then
        match = (str1 == str2)
      else
        call SplitDetails( str1, switch1, details1 )
        call SplitDetails( str2, switch2, details2 )
        match = (switch1 == switch2)
      endif
    end function matchem
  end subroutine GetUniqueStrings

  ! -------------------------------------------------  inserthashlement  -----
  subroutine INSERTHASHELEMENT ( NAME, VALUE, KEYS, VALUES, INSEPARATOR )
    ! Dummy args
    character(len=*), intent(in)           :: NAME
    character(len=*), intent(in)           :: VALUE
    character(len=*), intent(inout)        :: KEYS
    character(len=*), intent(inout)        :: VALUES
    character(len=1), intent(in), optional :: INSEPARATOR
    ! Local variables
    character(len=64) :: cvalue
    character (len=16)                            :: keyString
    character (len=8)                             :: nCh
    character (len=1)                             :: separator
    ! Executable
    separator = ','
    if ( present(inseparator) ) separator = inseparator
    ! 1st--is name an array-valued hash key?
    keyString = trim(name) // 'n'
    call GetHashElement( keys, values, keyString, nCh, &
      & countEmpty, inseparator=separator )
    if ( nCh == separator ) then
      ! No, it's just a scalar
      call PutHashElement ( keys, values, &
      & trim(name), value, countEmpty=countEmpty, inseparator=separator )
    else
      ! Yes, it's an array, so we must put it in two places:
      ! "name(cvalue)" and "name(n)" where
      ! cvalue is the actual value of "count"
      ! and "name(n)" is literally that (i.e., don't evaluate "n")
      call GetHashElement( keys, values, &
        & 'count', cvalue, countEmpty=countEmpty, inseparator=separator )
      keyString = trim(name) // '(n)'
      ! if ( DEEBUG ) then
      !  call outputnamedValue( 'keyString', trim(keyString) )
      !  call outputnamedValue( 'value', trim(value) )
      ! endif
      call PutHashElement ( keys, values, &
        & trim(keyString), value, countEmpty=countEmpty, inseparator=separator )
      keyString = trim(name) // '(' // trim(adjustl(cvalue)) // ')'
      ! if ( DEEBUG ) then
      !   call outputnamedValue( 'keyString', trim(keyString) )
      !   call outputnamedValue( 'value', trim(value) )
      ! endif
      call PutHashElement ( keys, values, &
        & trim(keyString), value, countEmpty=countEmpty, inseparator=separator )
    endif

  end subroutine INSERTHASHELEMENT

  ! -------------------------------------------------  Intersection  -----
  function Intersection ( str1, str2, options ) result ( outstr )
    ! return intersection of 2 string lists, blank means empty set
    ! E.g., given str1 = 'a,b,c' and str2 = 'd,e,f,c,a'
    ! returns 'a,c'
    ! options
    ! '-w' combined with wildcard '*' lets 'a*' match any string 'a..'
    !--------Argument--------!
    character (len=*), intent(in)                 :: str1
    character (len=*), intent(in)                 :: str2
    character (len=*), optional, intent(in)       :: options
    character (len=len(str1)+len(str2)+1)         :: outstr

    !----------Local vars----------!
    logical, parameter :: countEmpty = .true.
    logical, parameter :: ignoreLeadingSpaces = .true.
    character(len=len(str1)) :: elem
    character(len=len(str1)) :: uniq1
    character(len=len(str2)) :: uniq2
    integer :: n1, n2
    integer :: i
    !----------executable part----------!

    if ( str1 == ' ' .or. str2 == ' ' ) then
      outstr = ' '
      return
    endif
    
    call GetUniqueList( str1, uniq1, n1, &
      & ignoreLeadingSpaces=ignoreLeadingSpaces, options='-e' )
    call GetUniqueList( str2, uniq2, n2, &
      & ignoreLeadingSpaces=ignoreLeadingSpaces, options='-e' )
    outstr = ' '
    do i=1, n1
      call GetStringElement( uniq1, elem, i, countEmpty )
      if ( IsInList( uniq2, trim(elem), options) ) &
        & outstr = CatLists( outstr, elem )
    enddo

  end function Intersection

  ! ---------------------------------------------  IsInList  -----

  ! Is string in the stringList? options may expand criteria
  ! See notes above about options
  ! options
  ! Warning: wildcard may be found in either stringList or string
  ! which may or may not be what you intended
  ! E.g.
  ! stringList = 'abcd,a*,bcd,cd,d' and string='acd' with options = '-w'
  !              returns TRUE because 'a*' matches 'acd'
  ! Special cases:
  ! string == ' '      => always FALSE
  ! stringList == ' '  => always FALSE
  function IsInList( stringList, string, options ) result(itIs)
    ! Dummy arguments
    character (len=*), optional, intent(in)       :: stringlist
    character (len=*), intent(in)                 :: string
    character (len=*), optional, intent(in)       :: options
    logical                                       :: itIs
    ! Internal variables
    logical, parameter :: countEmpty = .true.
    integer :: i
    integer :: n
    character(len=MAXELEMENTLENGTH) :: element
    ! Executable
    itIs = .false.
    if ( .not. present(stringList) ) return
    if ( len_trim(stringList) < 1 .or. len_trim(string) < 1 ) return
    n = NumStringElements( stringList, countEmpty )
    do i=1, n
      call GetStringElement( stringList, element, i, countEmpty )
      if ( element == ' ' ) cycle
      itIs = streq( trim(element), trim(string), options )
      if ( itIs ) return
    enddo
  end function IsInList

  ! ---------------------------------------------  List2Array  -----

  ! This subroutine takes a (usually) comma-separated string list, interprets it
  ! as a list of individual elements and returns an equivalent array of
  ! sub-strings in which the n'th element is the n'th element

  ! if countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! if TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the separator may supplied, in case it isn't comma
  ! if the optional arg ignoreLeadingSpaces is TRUE, "a, b, c" is
  ! treated like "a,b,c"; otherwise the leading spaces are retained

  subroutine List2Array( inList, outArray, countEmpty, inseparator, &
   & IgnoreLeadingSpaces )
    ! Dummy arguments
    character (len=*), intent(in)                 :: inList
    character (len=*), dimension(:), intent(out)  :: outArray
    logical, intent(in)                           :: countEmpty
    character (len=*), optional, intent(in)       :: inseparator
    logical, optional, intent(in)                 :: IgnoreLeadingSpaces

    ! Local variables
    integer :: elem, nElems

    logical                         :: myIgnoreLeadingSpaces
    ! Executable code

    if(present(IgnoreLeadingSpaces)) then
      myIgnoreLeadingSpaces = IgnoreLeadingSpaces
    else
      myIgnoreLeadingSpaces = .false.
    endif

    if ( size(outArray) <= 0 ) return
    outArray = BLANK
    elem = 1
    nElems = NumStringElements(inList, countEmpty, inseparator)
    if ( nElems <= 0 ) return
    do
      call GetStringElement(inList, outArray(elem), elem, countEmpty, inseparator)
      if ( myIgnoreLeadingSpaces ) outArray(elem) = adjustl(outArray(elem))
      elem = elem + 1
      if ( elem > min(nElems, size(outArray)) ) return
    enddo

  end subroutine List2Array

  ! ---------------------------------------------  listMatches  -----

  ! Return list of matches for string in the stringList. 
  ! options may expand criteria
  ! See notes above about options
  ! options
  ! Warning: wildcard may be found in either stringList or string
  ! which may or may not be what you intended
  ! E.g.
  ! stringList = 'abcd,a*,bcd,cd,d' and string='acd' with options = '-w'
  !              returns 'a*' because 'a*' matches 'acd'
  ! Special cases:
  ! string == ' '      => always ' '
  ! stringList == ' '  => always ' '
  function listMatches( stringList, string, options ) result(matches)
    ! Dummy arguments
    character (len=*), intent(in)                 :: stringlist
    character (len=*), intent(in)                 :: string
    character (len=*), optional, intent(in)       :: options
    character (len=len(stringList))               :: matches
    ! Internal variables
    ! logical, parameter :: countEmpty = .true.
    integer :: i
    logical :: itMatches
    integer :: n
    character(len=max(len(stringList), len(string))) :: element
    ! Executable
    call prepOptions( options )
    matches = ' '
    if ( len_trim(stringList) < 1 .or. len_trim(string) < 1 ) return
    n = NumStringElements( stringList, countEmpty )
    do i=1, n
      call GetStringElement( stringList, element, i, countEmpty )
      if ( element == ' ' ) cycle
      itMatches = streq( trim(element), trim(string), options )
      if ( itMatches ) then
        matches = CatLists( matches, element )
      endif
    enddo
  end function listMatches

  ! -------------------------------------------------  LoopOverFormula  -----
  subroutine LoopOverFormula ( formula, arg, values, results )
    ! Looping while it evaluates a formula, plugging in the nth value for
    ! each occurrence of the arg appearing as '${arg}'
    ! E.g., if arg is "phase", formula is 
    ! "x${phase}: vector, template=state"
    ! and the values of arg are (/"InitPtan ", "FinalPtan"/) 
    ! then the results array will be
    ! (/
    ! "xInitPtan: vector, template=state", 
    ! "xFinalPtan: vector, template=state"
    ! /)
    !--------Argument--------!
    character (len=*), intent(in)                :: formula
    character (len=*), intent(in)                :: arg
    character (len=*), dimension(:), intent(in)  :: values
    character (len=*), dimension(:), intent(out) :: results

    !----------Local vars----------!
    integer :: i, n
    character (len=len(formula)+len(values)) :: tmpstr
    character (len=len(values))              :: value
    character(len=len_trim(arg)+8)           :: variable
    !----------executable part----------!
    n = min( size(values), size(results) )
    results = formula
    ! Check whether we have any work to do
    if ( n < 1 .or. index( formula, '${' ) < 1 ) return
    variable = '${' // trim(adjustl(arg)) // '}'
    do i=1, n
      value = values(i)
      tmpstr = results(i)
      call ReplaceSubString( tmpStr, results(i), trim(variable), trim(value), &
        & which='all', no_trim=.true. )
    enddo
  end subroutine LoopOverFormula

  ! .............................................  nCharsinFormat  .....
  function nCharsinFormat ( Format ) result(nplusm)
    ! Utility to calculate how many characters in a format spec:         
    ! [n{xX}][,]{DEFGdefg}m.b                                             
    ! where n, m, and b are digits (we care only about n and m)           
    ! return (n+m)
    ! Tested for specs: sci. format esm.b and eng. format enm.b
    ! Also for min. width spec: 'f0.b' it will silently return 0
    ! (It's up to you to handle that correctly)
    ! Args                                                                
    character(len=*), intent(in) ::  Format                               
    integer :: nplusm                                                     
    ! Local variables                                                     
    character(len=20) :: kChar, myFormat                                  
    integer :: n, m
    ! Executable                                                          
    nplusm = 0                                                            
    kChar=lowerCase(Format)
    call ReplaceSubString(kChar, myFormat, 'es', 'f')                   
    call ReplaceSubString(myFormat, kChar, 'en', 'f')                   
    call ReplaceSubString(kChar, myFormat, 'g', 'f')                   
    call ReplaceSubString(myFormat, kChar, 'e', 'f')                   
    call ReplaceSubString(kChar, myFormat, 'd', 'f')                   
    call ExtractSubString(TRIM(myFormat), kChar, 'f', '.')             
    if ( kChar == '0' ) return ! Special case of e.g. 'f0.3'
    read (kChar, '(i2)') m                                                
    if (m < 1) call PrintMessage ( MLSMSG_Error, ModuleName, &              
      & 'Bad conversion to m in OUTPUT_xxxLE (format not "{defg}"' )      
    if ( index(TRIM(myFormat), 'x' ) == 0 ) then                          
      n = 0                                                               
    else                                                                  
      call ExtractSubString(TRIM(myFormat), kChar, '(', 'x')           
      read (kChar, '(i2)') n                                              
      if (n < 1) then                                                     
        print *, trim(kChar)                                              
        print *, trim(myFormat)                                           
        call PrintMessage ( MLSMSG_Error, ModuleName, &                     
          & 'Bad conversion to n in OUTPUT_xxxLE (format not "{defg}"' )  
      end if                                                              
    end if                                                                 
    nplusm = n + m                                                        
  end function nCharsinFormat

  ! ---------------------------------------------  NumStringElements  -----

  ! This function takes a (usually) comma-separated string list, interprets it
  ! as a list of individual elements and returns the
  ! number of elements
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-separated lists
  !
  ! if countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE  

  ! As an optional arg the separator may supplied, in case it isn't comma

  ! See also GetStringElement

  function NumStringElements(inList, countEmpty, &
   & inseparator, Longestlen) RESULT (nElements)
    ! Dummy arguments
    character (len=*), intent(in)             :: inList
    logical, intent(in)                       :: countEmpty
    integer                                   :: nElements
    character (len=*), optional, intent(in)       :: inseparator
    integer, optional, intent(out)            :: Longestlen  ! LENGTH of longest

    ! Local variables
    integer :: i, sinceLastseparated           ! Loop counters
    logical :: lastWasNotseparated

    character (len=1)               :: separator
    ! Executable code

    if(present(inseparator)) then
      separator = inseparator
    else
      separator = COMMA
    endif

   ! Count the number of separators
   if ( present(Longestlen) ) &
     & Longestlen =0
   ! nElements-1 = number of separators
   if(len_TRIM(inList) <= 0) then
     nElements=0
      if ( present(Longestlen) ) Longestlen = 0
      RETURN
   endif

   lastWasNotseparated = .FALSE.
   nElements = 1
   sinceLastseparated = 0
   do i=1, len_TRIM(inList)
     if(inList(i:i) == separator) then
       if(countEmpty .OR. lastWasNotseparated) then
         nElements = nElements+1
            if ( present(Longestlen) ) &
             & Longestlen = max(Longestlen, sinceLastseparated)
       endif
       lastWasNotseparated = .FALSE.
       sinceLastseparated = 0
     else
       lastWasNotseparated = .TRUE.
       sinceLastseparated = sinceLastseparated + 1
     endif
   enddo
   if ( present(Longestlen) ) &
     & Longestlen = max(Longestlen, sinceLastseparated)

  end function NumStringElements

  ! ---------------------------------------------  optionDetail  -----

  ! This function takes a string, interprets it
  ! as a list of individual one-character options, and a test option
  ! It returns
  !  'yes' if the option is present
  !  'no'  if the option is absent
  !  'arg' if the option is present and completed by (takes as argument) 'arg'
  !    
  ! single_option is a one-character option
  ! multi_option is a multi-character alternate form
  ! E.g., single_option might be 'a' while multi_option might be 'answer'
  ! so the option could be set in either form '-a' or '--answer'
  ! if patterns is present, it determines whether options must be preceded
  ! by an '-' and how args are to be denoted
  ! Recognized values of pattern are
  ! 0: '-ab[arg] --xyz=arg' means  '(default)'
  !    may catenate single-char single_options; any arg is surrounded by "[]"
  !    multiple-char multi_option preceded by '--'; any arg set off by "="
  ! 1: '-a -b arg --xyz=arg' means 
  !    each single_option preceded by '-'; any arg is set off by a space
  !    multiple-char multi_option preceded by '--'; any arg set off by "="
  ! 2: '-a -b arg -xyz arg' means 
  !    each single_option preceded by '-'; any arg is set off by a space
  !    multiple-char multi_option preceded by '-'; any arg set off by a space
  ! 3: '-a -b arg --xyz arg' means 
  !    each single_option preceded by '-'; any arg is set off by a space
  !    multiple-char multi_option preceded by '--'; any arg set off by a space
  
  ! 4: '-a -barg --xyz arg' means 
  !    each single_option preceded by '-'; followed immediately by any arg
  !    multiple-char multi_option preceded by '--'; any arg set off by a space
  
  ! As an example, say the list of options is
  ! "-ab[arg1]c[arg2]d"
  ! and the test option is "b"
  ! The returned value would be 'arg1'
  ! if the test option were "a" the returned value would be 'yes'
  ! if the test option were "c" the returned value would be 'arg2'
  ! if the test option were "g" the returned value would be 'no'
  ! (because g doesn't count except outside the '[]' chars)
  
  ! The behavior may be modified by pattern and delims args
  ! For which see comment above
  
  ! Notes:
  ! (1) if the string list is absent then the test option is automatically absent
  ! (2) if the string list is "*" then
  ! the test option is automatically present (would you like to override that?)
  ! (3) Why don't you let the '[]' pair that set off args be
  ! overridden, say by other optional args?
  
  function optionDetail( inList, single_option, multi_option, &
    & pattern, delims ) RESULT (detail)
    ! Dummy arguments
    character (len=*), optional, intent(in)   :: inlist
    character (len=1), optional, intent(in)   :: single_option
    character (len=*), optional, intent(in)   :: multi_option
    integer, optional, intent(in)             :: pattern
    character (len=1), dimension(2), optional, &
     & intent(in)                             :: delims
    character (len=MAXELEMENTLENGTH)               :: detail
    ! Local variables
    integer :: bloc
    logical, parameter :: COUNTEMPTY = .true.
    character :: cquotes, quotes
    integer :: k
    character (len=MAXELEMENTLENGTH)           :: element
    character (len=MAXELEMENTLENGTH)           :: listBloc ! space-separated
    logical :: multi
    integer :: myPattern
    character(len=16) :: test_multi
    character :: test_option
    character(len=*), dimension(0:4), parameter :: multi_prefix = &
      & (/ '--', '--', '- ', '--' , '--' /)

    ! Executable code
    detail = 'no'
    if ( .not. present(inList) ) return
    test_option = char(0) ! NULL
    test_multi = char(0) ! NULL
    if ( present(single_option) ) test_option = single_option
    if ( test_option == ' ' ) test_option = char(0) ! NULL
    if ( present(multi_option) ) test_multi = multi_option
    if ( test_multi == ' ' ) test_multi = char(0) ! NULL
    myPattern = 0
    if ( present(pattern) ) then
      if ( any(Pattern == (/0, 1, 2, 3, 4 /)) &  ! accept legal values only
        & ) myPattern = pattern
    endif
    if ( adjustl(inList) == '*' ) then
      detail = 'yes'
      return
    endif
    quotes = '['
    cquotes = ']'
    if ( present(delims) ) then
      quotes = delims(1)
      cquotes = delims(2)
    endif
    if ( .not. present(single_option) .and. .not. present(multi_option) ) return
    select case (myPattern)
    case ( 0 )
      ! '-ab[arg] --xyz=arg
      if ( index( inList, test_option ) < 1 ) then
        if ( .not. present(multi_option) ) return
      endif

      ! OK, test_option or its alt may be present, but where? 
      ! does it have an arg?

      do bloc = 1, NumStringElements( inList, countEmpty, inseparator=' ' )
        listBloc = StringElement( inList, bloc, countEmpty, inseparator=' ' )
        ! does this block begin with one "-" or two?
        multi = ( index( listBloc, trim(multi_prefix(myPattern)) ) > 0 )
        if ( multi ) then
          if ( .not. present(multi_option) ) cycle
          k = index(listBloc, trim(multi_prefix(myPattern)) // &
            & trim(test_multi) // '=' )
          if ( k > 0 ) then
            k = index(listBloc, '=' )
            detail = listBloc(k+1:)
            return
          elseif ( index(listBloc, trim(multi_prefix(myPattern)) // &
            & trim(test_multi) // ' ' ) > 0 ) then
            detail = 'yes'
            return
          endif
        else
          ! 1st--rid ourselves of everything bracketed by '[]'
          ! element = unquote( listBloc, quotes='[', cquotes=']', &
          element = unquote( listBloc, quotes=quotes, cquotes=cquotes, &
            & options='-r' )
          ! print *, 'After unquote: ', trim(element)
          if ( index(element, test_option) < 1 ) cycle
          call extractSubstring( listBloc, element, &
            & test_option // quotes, cquotes )
          ! print *, 'After extracting: ', trim(element)
          if ( len_trim(element) > 0 ) then
            detail = element
          else
            detail = 'yes'
          endif
          return
        endif
      enddo

    case ( 1 )
      ! '-a -b arg --xyz=arg'
      if ( index( inList, test_option ) < 1 ) then
        if ( .not. present(multi_option) ) return
      endif

      ! OK, test_option or its alt may be present, but where? 
      ! does it have an arg?

      do bloc = 1, NumStringElements( inList, countEmpty, inseparator=' ' )
        listBloc = StringElement( inList, bloc, countEmpty, inseparator=' ' )
        ! does this block begin with one "-" or two?
        multi = ( index( listBloc, trim(multi_prefix(myPattern)) ) > 0 )
        if ( multi ) then
          if ( .not. present(multi_option) ) cycle
          k = index(listBloc, trim(multi_prefix(myPattern)) // &
            & trim(test_multi) // '=' )
          if ( k > 0 ) then
            k = index(listBloc, '=' )
            detail = listBloc(k+1:)
            return
          elseif ( index(listBloc, trim(multi_prefix(myPattern)) // &
            & trim(test_multi) // ' ' ) > 0 ) then
            detail = 'yes'
            return
          endif
        else
          if ( index(listBloc, '-' // test_option) > 0 ) then
            ! OK, we've got the option all right; but is the next block an arg?
            detail = 'yes'
            if ( bloc == NumStringElements( inList, countEmpty, inseparator=' ' ) ) &
              & return
            element = StringElement( inList, bloc+1, countEmpty, inseparator=' ' )
            if ( index(adjustl(element), '-') == 1 ) then
              return
            else
              detail = element
              return
            endif
          endif
        endif
      enddo

    case ( 2, 3 )
      ! '-a -b arg -xyz arg' or
      ! '-a -b arg --xyz arg'
      if ( index( inList, test_option ) < 1 ) then
        if ( .not. present(multi_option) ) return
      endif

      ! OK, test_option or its alt may be present, but where? 
      ! does it have an arg?

      do bloc = 1, NumStringElements( inList, countEmpty, inseparator=' ' )
        listBloc = StringElement( inList, bloc, countEmpty, inseparator=' ' )
        ! does this block begin with one "-" or not?
        if ( index(adjustl(listbloc), '-') == 1 ) then
          if ( adjustl(listBloc) == '-' // test_option ) then
            detail='yes' ! keep going--next we'll check for an arg
          elseif( .not. present(multi_option) ) then
            cycle
          elseif ( adjustl(listBloc) == trim(multi_prefix(myPattern)) // &
            & test_multi ) then
            detail='yes' ! keep going--next we'll check for an arg
          else
            cycle
          endif
          ! Now check if next bloc is an arg
          if ( bloc == NumStringElements( inList, countEmpty, inseparator=' ' ) ) &
            & return
          element = StringElement( inList, bloc+1, countEmpty, inseparator=' ' )
          if ( index(adjustl(element), '-') == 1 ) then
            return
          else
            detail = element
            return
          endif
        endif
      enddo
    case ( 4 )
      ! '-a -barg --xyz arg'
      if ( index( inList, '-' // test_option ) < 1 ) then
        if ( .not. present(multi_option) ) return
      endif

      ! OK, test_option or its alt may be present, but where? 
      ! does it have an arg?

      do bloc = 1, NumStringElements( inList, countEmpty, inseparator=' ' )
        listBloc = StringElement( inList, bloc, countEmpty, inseparator=' ' )
        ! does this block begin with one "-" or two?
        multi = ( index( listBloc, trim(multi_prefix(myPattern)) ) > 0 )
        if ( multi ) then
          if ( .not. present(multi_option) ) cycle
          if ( index( listBloc, trim(multi_prefix(myPattern)) // &
            & trim(multi_option) ) < 1 ) cycle
          detail = StringElement( inList, bloc+1, countEmpty, inseparator=' ' )
          if ( detail == ' ' .or. detail(1:1) == '-' ) detail = 'yes'
        else
          if ( index(listBloc, '-' // test_option) > 0 ) then
            ! OK, we've got the option all right; but is it followed by an arg?
            detail = 'yes'
            if ( listBloc(3:3) == ' ' ) cycle
            detail = listBloc(3:)
          endif
        endif
      enddo
    case default
    end select
  end function optionDetail

  ! ---------------------------------------------  ParseOptions  -----

  ! Parse a commandline, checking for the presence of single-character
  ! options or multi-character options, returning in opts_out for each
  !  'yes' if the option is present
  !  'no'  if the option is absent
  !  'arg' if the option is present and completed by its needed 'arg'
  !  ' '   if the option is present but missing its needed 'arg'
  !    
  ! single_options are all the one-character options
  ! multi_options are multi-character alternate forms
  ! whether it needs an arg is determined by needs_arg
  
  ! cmdargs, if present, will return any remaining commandline args

  ! The behavior may be modified by pattern and delims args
  ! For which see comment above
  
  ! Notes:
  ! (1) See also optionDetail, switchDetail
  
  subroutine ParseOptions( cmdline, opts_out, pattern, single_options, &
    & multi_options, needs_arg, delims, cmdargs )
    ! Dummy arguments
    character (len=*), intent(in)                  :: cmdline
    character (len=*), dimension(:), intent(out)   :: opts_out
    character (len=1), dimension(:), intent(in)    :: single_options
    integer, intent(in)                            :: pattern
    character (len=*),  dimension(:), intent(in)   :: multi_options
    logical,  dimension(:), intent(in)             :: needs_arg
    character (len=1), dimension(2), optional, &
      & intent(in)                                 :: delims
    character (len=*), dimension(:), optional, intent(out)&
      &                                            :: cmdargs

    ! Local variables
    logical, parameter :: COUNTEMPTY = .true.
    integer :: i, j, k
    character(len=16)  :: lastOption
    logical            :: lastOptionNeededArg
    integer :: nopts
    ! Begin executable
    if ( size(opts_out) > 0 ) opts_out = 'no'
    if ( present(cmdargs) ) cmdargs = ' '
    if ( size(opts_out) < 1 ) return
    if ( len_trim(cmdline) < 1 ) return
    nopts = min( size(single_options), size(multi_options) )
    nopts = min( nopts, size(opts_out) )
    do i = 1, nopts
      opts_out(i) = optionDetail( cmdline, single_options(i), &
        & multi_options(i), pattern, delims )
      ! Now check for needed args
      if ( needs_arg(i) ) then
        if ( opts_out(i) == 'yes' ) opts_out(i) = ' '
      else
        if ( .not. IsInList('yes,no', trim(opts_out(i))) ) opts_out(i) = 'yes'
      endif
    enddo
    if ( .not. present(cmdargs) ) return
    ! Must find last option
    k = index( cmdline, '-', back=.true. )
    if ( k < 1 ) then
      ! No options on cmdline--everything is a cmdarg
      do i = 1, NumStringElements( cmdline, countEmpty, inseparator=' ' )
        cmdargs(i) = StringElement( cmdline, i, countEmpty, inseparator=' ' )
      enddo
      return
    endif
    ! Now must find which option that was
    lastOption = StringElement( cmdline(k+1:), 1, countEmpty, inseparator=' ' )
    i = FindFirst( single_options, trim(lastOption) )
    if ( i < 1 ) i = FindFirst( multi_options, trim(lastOption) )
    if ( i < 1 ) then
      ! Whoa, option not recognized, but we'll just carry on nonetheless
      lastOptionNeededArg = .false.
    else
      lastOptionNeededArg = needs_arg(i)
    endif
    i = 0
    do j=1, NumStringElements( cmdline(k+1:), countEmpty, inseparator=' ' )
      if ( j == 1 ) cycle ! skip last option
      if ( j == 2 .and. lastOptionNeededArg ) cycle ! skip last option's arg
      if ( i >= size(cmdargs) ) cycle
      i = i + 1
      cmdargs(i) = StringElement( cmdline(k+1:), j, countEmpty, inseparator=' ' )
    enddo
  end subroutine ParseOptions

  ! ---------------------------------------------  PutHashElement  -----

  ! This family of subroutines interprets two arguments as
  ! a set of {key = value} pairs
  ! Given a (possibly new) key and value, insert or replace the value
  ! This subroutine takes two (usually) comma-separated string lists, interprets it
  ! each as a list of elements, treating the first as keys and the second as
  ! a hash table, associative array or dictionary
  ! It replaces with elem the sub-string from the hash table corresponding to the key
  ! if the key is not found in the array of keys, it adds a new key

  ! if countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! if TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the separator may supplied, in case it isn't comma
  ! Another optional arg, part_match, returns a match for the 
  ! first hash element merely found in the key; e.g.
  ! 'won, to, tree' and key 'protocol.dat' matches 'to'

  ! Basic premise: Find the element number corresponding to the key
  ! if found remove that element from both key and hash list
  ! then add new key and hash to lists

  ! Someday you may wish to define a StringHash_T made up of the two
  ! strings
  

  subroutine PutHashElement_int( keys, values, key, value, &
  & countEmpty, inseparator, part_match )
    ! Dummy arguments
    character (len=*), intent(inout)          :: keys
    integer, dimension(:), intent(inout)      :: values
    character (len=*), intent(in)             :: key
    integer, intent(in)                       :: value
    logical, intent(in)                       :: countEmpty
    character (len=*), optional, intent(in)   :: inseparator
    logical, optional, intent(in)             :: part_match

    ! Local variables
    integer                                    :: N
    integer                                    :: num

    ! Executable code

    num = StringElementNum(keys, key, countEmpty, inseparator, part_match)
    if( num > 0 .and. num <= size(values) ) then
      values(num) = value
    elseif( num > size(values) ) then
      ! Can't handle arrays this big
    else
      ! key not found :: must add to keys, values
      N = NumStringElements( keys, countEmpty, inseparator )
      keys = CatLists( keys, key, inseparator )
      values(N+1) = value
    endif

  end subroutine PutHashElement_int

  subroutine PutHashElement_log( keys, values, key, value, &
  & countEmpty, inseparator, part_match )
    ! Dummy arguments
    character (len=*), intent(inout)          :: keys
    logical, dimension(:), intent(inout)      :: values
    character (len=*), intent(in)             :: key
    logical, intent(in)                       :: value
    logical, intent(in)                       :: countEmpty
    character (len=*), optional, intent(in)   :: inseparator
    logical, optional, intent(in)             :: part_match

    ! Local variables
    integer                                    :: N
    integer                                    :: num

    ! Executable code
    num = StringElementNum( keys, key, countEmpty, inseparator, part_match )
    if( num > 0 .and. num <= size(values) ) then
      values(num) = value
    elseif( num > size(values) ) then
      ! Can't handle arrays this big
    else
      ! key not found :: must add to keys, values
      N = NumStringElements( keys, countEmpty, inseparator )
      keys = CatLists( keys, key, inseparator )
      values(N+1) = value
    endif

  end subroutine PutHashElement_log

  subroutine PutHashElement_str( keyList, hashList, key, elem, &
  & countEmpty, inseparator, part_match )
    ! Dummy arguments
    character (len=*), intent(inout)          :: keyList
    character (len=*), intent(inout)          :: hashList
    character (len=*), intent(in)             :: key
    character (len=*), intent(in)             :: elem
    logical, intent(in)                       :: countEmpty
    character (len=*), optional, intent(in)   :: inseparator
    logical, optional, intent(in)             :: part_match

    ! Local variables
    integer                                    :: num
    character(len=len(keyList)+len(key)+1)     :: keys
    character(len=len(hashList)+len(elem)+1)   :: hash

    ! Executable code
    num = StringElementNum( keyList, key, countEmpty, inseparator, part_match )
    if( num > 0) then
      call RemoveNumFromList( keyList, keys, num, inseparator )
      call RemoveNumFromList( hashList, hash, num, inseparator )
    else
      keys = keyList
      hash = hashList
    endif
    keyList = CatLists( keys, key, inseparator )
    hashList = CatLists( hash, elem, inseparator )

  end subroutine PutHashElement_str

  subroutine PutHashElement_strarray( KEYLIST, HASHLIST, KEY, ARRAY, &
  & COUNTEMPTY, INSEPARATOR, PART_MATCH )
    ! We insert an array of values onto a hash
    ! storing the array like this
    ! array "name" contains "value_1", "value_2" .. "nalue_n"
    !    key        value
    !    ---        -----
    ! "namen"        "n" (number of elements)
    ! "name(1)"   "value_1"
    ! "name(2)"   "value_2"
    !    .    .    .
    ! "name(n)"   "value_n"
    ! Dummy arguments
    character (len=*), intent(inout)              :: KEYLIST
    character (len=*), intent(inout)              :: HASHLIST
    character (len=*), intent(in)                 :: KEY
    character (len=*), dimension(:), intent(in)   :: ARRAY
    logical, intent(in)                           :: COUNTEMPTY
    character (len=*), optional, intent(in)       :: INSEPARATOR
    logical, optional, intent(in)                 :: PART_MATCH

    ! Local variables
    integer                                       :: j
    character (len=16)                            :: keyString
    integer                                       :: n
    character (len=8)                             :: nCh

    ! Executable code
    n = size(array)
    if ( n < 1 ) return
    call writeIntsToChars( n, nCh )
    keyString = trim(key) // 'n'
    call PutHashElement_str( keyList, hashList, keyString, nCh, &
      & countEmpty, inseparator, part_match )
    
    do j=1, n
      call writeIntsToChars( j, nCh )
      keyString = trim(key) // '(' // trim(adjustl(nCh)) // ')'
      call PutHashElement_str( keyList, hashList, keyString, array(j), &
        & countEmpty, inseparator, part_match )
    enddo

  end subroutine PutHashElement_strarray

  ! --------------------------------------------------  ReadIntsFromList  -----
  subroutine ReadIntsFromList ( inList, ints, error )
    ! Takes a list and reads it as an array of ints
    ! E.g., given '1 2 2 3 4 5'  returns (/ 1, 2, 2, 3, 4, 5 /)
    ! (Inverse of WriteIntsToList)
    !--------Argument--------!
    character (len=*), intent(in)      :: inList
    integer, dimension(:), intent(out) :: ints
    integer, optional, intent(out)     :: error
    ! Method:
    ! Use Fortran read
    integer :: status
    ints = -999
    status = 0
    if ( len_trim(inList) > 0 ) read(inList, *, iostat=status, err=100, end=100) ints
100   continue
    if ( present(error) ) error = status
  end subroutine ReadIntsFromList

  ! --------------------------------------------------  ReadNumsFromList  -----
  subroutine ReadDoubleArrayFromString ( inList, nums, separator, ignore, error )
    ! Takes a list and reads it as an array of floats
    ! E.g., given '6.32 0 9.05'  returns (/ 6.32, 0., 9.05 /)
    ! Ignores any non-numerical stuff found in ignore if you supply ignore
    !--------Argument--------!
    character (len=*), intent(in)                        :: inList
    double precision, dimension(:), intent(out)          :: nums
    integer, optional, intent(out)                       :: error
    character (len=*), optional, intent(in)              :: separator
    character (len=*), optional, intent(in)              :: ignore
    ! Method:
    ! (1) Turn list to an array
    ! (2) Read each array element separately into a float
    integer                                              :: nElems
    character(len=16), dimension(MAXSTRELEMENTLENGTH)    :: strArray
    !
    if ( present(error) ) error = 0
    call List2Array( inList, strArray, &
      & countEmpty=.false., inseparator=separator )
    nElems = FindLast ( strArray, ' ', reverse=.true. )
    if ( nElems < 1 ) then
      if ( present(error) ) error = -999
      return
    endif
    call readNumsFromChars( strArray(:nElems), nums(:nElems), ignore=ignore )
  end subroutine ReadDoubleArrayFromString

  subroutine ReadRealArrayFromString ( inList, nums, separator, ignore, error )
    ! Takes a list and reads it as an array of floats
    ! E.g., given '6.32 0 9.05'  returns (/ 6.32, 0., 9.05 /)
    ! Ignores any non-numerical stuff found in ignore if you supply ignore
    !--------Argument--------!
    character (len=*), intent(in)                        :: inList
    real, dimension(:), intent(out)                      :: nums
    integer, optional, intent(out)                       :: error
    character (len=*), optional, intent(in)              :: separator
    character (len=*), optional, intent(in)              :: ignore
    ! Method:
    ! (1) Turn list to an array
    ! (2) Read each array element separately into a float
    integer                                              :: nElems
    character(len=16), dimension(MAXSTRELEMENTLENGTH)    :: strArray
    !
    if ( present(error) ) error = 0
    call List2Array( inList, strArray, &
      & countEmpty=.false., inseparator=separator )
    nElems = FindLast ( strArray, ' ', reverse=.true. )
    if ( nElems < 1 ) then
      if ( present(error) ) error = -999
      return
    endif
    call readNumsFromChars( strArray(:nElems), nums(:nElems), ignore=ignore )
  end subroutine ReadRealArrayFromString

  ! --------------------------------------------------  RemoveElemFromList  -----
  subroutine RemoveElemFromList ( inList, outList, elem, inseparator, &
    & options )
    ! Takes a list and removes all occurrence(s) of elem
    ! E.g., given 'a,b,c,d,..,z' and asked to remove 'c' returns 'a,b,d,..z'
    !--------Argument--------!
    character (len=*), intent(in) :: inList
    character (len=*), intent(in) :: elem
    character (len=*), intent(out)                :: outList
    character (len=*), optional, intent(in)       :: inseparator
    character (len=*), optional, intent(in)       :: options
    ! Method:
    ! Prepend elem onto start of list, make it unique,
    ! then snip it back off
    !----------Local vars----------!
    logical :: myCountEmpty
    character(len=8) :: myOptions
    integer :: numUnique
    character(len=len(inList)+len(elem)+1) :: temp_list, unique_list
    character (len=1)               :: separator
    !----------Executable part----------!
    if(present(inseparator)) then
      separator = inseparator
    else
      separator = COMMA
    endif
    myOptions = ' '
    if ( present(options) ) myOptions = options
    myCountEmpty = index( myOptions, 'e' ) > 0  ! .true.

    outList = inList
    if (len_trim(elem) < 1 .or. len_trim(inList) < 1 &
      & .or. StringElementNum(inList, elem, myCountEmpty, &
    & inseparator=inseparator) < 1 ) RETURN
    temp_list = trim(elem) // separator // trim(inList)
    call GetUniqueList( temp_list, unique_list, numUnique, &
    & inseparator=inseparator, ignoreLeadingSpaces=.true., options=options )
    ! outList = unique_list(len(elem)+1:)
    ! The following is evidence of poor programming habits
    ! (As if any more evidence was needed)
    if ( unique_list(len_trim(elem)+1:len_trim(elem)+1) == separator ) then
      outList = unique_list(len_trim(elem)+2:)
    else
      outList = unique_list(len_trim(elem)+1:)
    endif
  end subroutine RemoveElemFromList

  ! ---------------------------------------------  RemoveHashArray  -----
  subroutine RemoveHashArray( KEYLIST, HASHLIST, KEY, &
  & COUNTEMPTY, INSEPARATOR, PART_MATCH )
    ! We remove an array of values from a hash
    ! storing the array like this
    ! array "name" contains "value_1", "value_2" .. "nalue_n"
    !    key        value
    !    ---        -----
    ! "namen"        "n" (number of elements)
    ! "name(1)"   "value_1"
    ! "name(2)"   "value_2"
    !    .    .    .
    ! "name(n)"   "value_n"
    ! Dummy arguments
    character (len=*), intent(inout)              :: KEYLIST
    character (len=*), intent(inout)              :: HASHLIST
    character (len=*), intent(in)                 :: KEY
    logical, intent(in)                           :: COUNTEMPTY
    character (len=*), optional, intent(in)       :: INSEPARATOR
    logical, optional, intent(in)                 :: PART_MATCH

    ! Local variables
    integer                                       :: j
    character (len=16)                            :: keyString
    integer                                       :: n
    character (len=8)                             :: nCh
    character (len=1)                             :: separator

    ! Executable code

    if(present(inseparator)) then
      separator = inseparator
    else
      separator = COMMA
    endif
    
    keyString = trim(key) // 'n'
    call GetHashElement_str( keyList, hashList, keyString, nCh, &
      & countEmpty, inseparator, part_match )
    if ( nCh == separator ) return
    call readIntsFromChars( nCh, n )
    call RemoveHashElement_str( keyList, hashList, keyString, &
      & countEmpty, inseparator, part_match )
    do j=1, n
      call writeIntsToChars( j, nCh )
      keyString = trim(key) // '(' // trim(adjustl(nCh)) // ')'
      call RemoveHashElement_str( keyList, hashList, keyString, &
        & countEmpty, inseparator, part_match )
    enddo

  end subroutine RemoveHashArray

  ! ---------------------------------------------  RemoveHashElement  -----

  ! This subroutine removes the key and corresponding value from a hash
  ! if key is not found among the keyList, it does nothing
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-separated lists

  ! if countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! if TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the separator may supplied, in case it isn't comma
  ! See also SplitWords

  subroutine RemoveHashElement_str( keyList, hashList, key, &
  & countEmpty, inseparator, part_match )
    ! Dummy arguments
    character (len=*), intent(inout)          :: keyList
    character (len=*), intent(inout)          :: hashList
    character (len=*), intent(in)             :: key
    logical, intent(in)                       :: countEmpty
    character (len=*), optional, intent(in)   :: inseparator
    logical, optional, intent(in)             :: part_match

    ! Local variables
    integer                                    :: num
    character(len=len(keyList)+1)              :: keys
    character(len=len(hashList)+1)             :: hash

    ! Executable code
    num = StringElementNum( keyList, key, countEmpty, inseparator, part_match )
    if( num > 0) then
      call RemoveNumFromList( keyList, keys, num, inseparator )
      call RemoveNumFromList( hashList, hash, num, inseparator )
      keyList = keys
      hashList = hash
    endif

  end subroutine RemoveHashElement_str

  ! ------------------------------------------------  RemoveListFromList  -----
  subroutine RemoveListFromList ( inList, outList, exclude, &
    & inseparator, options )
    ! Takes one list and removes from it all occurrence(s) 
    ! of each elem in another list called "exclude"
    ! E.g., given 'a,b,c,d,..,z' and asked to remove 'c,a' returns 'b,d,..z'
    !--------Argument--------!
    character (len=*), intent(in) :: inList
    character (len=*), intent(in) :: exclude ! What to exclude
    character (len=*), intent(out)                :: outList
    character (len=*), optional, intent(in)       :: inseparator
    character (len=*), optional, intent(in)       :: options
    ! Method:
    ! Repeatedly call RemoveElemFromList for each elem of exclude
    !----------Local vars----------!
    integer :: i
    character(len=max(len(inList), len(exclude)) + 1) :: elem
    logical :: myCountEmpty
    character(len=8) :: myOptions
    integer :: numElems
    character(len=len(inList)+1) :: temp_list
    character(len=len(inList)+1) :: temp_list2
    logical :: verbose
    !----------Executable part----------!
    myOptions = ' '
    if ( present(options) ) myOptions = options
    myCountEmpty = index( myOptions, 'e' ) > 0  ! .true.
    verbose = index( myOptions, 'v' ) > 0  ! .true.
    ! if ( present(countEmpty) ) myCountEmpty = countEmpty
    outList = inList
    if ( len_trim(exclude) < 1 .or. len_trim(inList) < 1 ) return
    numElems = NumStringElements( exclude, myCountEmpty, inseparator )
    if ( verbose ) print *, 'numElems in exclude ', numElems
    if ( numElems < 1 ) return
    temp_list = inList
    do i=1, numElems
      call GetStringElement( exclude, elem, i, myCountEmpty, inseparator )
      call RemoveElemFromList( temp_list, temp_list2, elem, &
        & inseparator, options )
      if ( verbose ) then
        print *, 'After excluding ' // trim(elem)
        print *, 'list was ' // trim(temp_list)
        print *, 'now ' // trim(temp_list2)
      endif
      temp_list = temp_list2
    enddo
    outList = temp_list
  end subroutine RemoveListFromList

  ! --------------------------------------------------  RemoveNumFromList  -----
  subroutine RemoveNumFromList ( inList, outList, nElement, inseparator, &
    & options )
    ! Removes a numbered element from a list
    ! E.g., given 'a,b,c,d,..,z' and asked to remove number 3 returns 'a,b,d,..z'
    !--------Argument--------!
    character (len=*), intent(in) :: inlist
    integer          , intent(in) :: nElement
    character (len=*), intent(out)                :: outList
    character (len=*), optional, intent(in)       :: inseparator
    character (len=*), optional, intent(in)       :: options
    ! Method:
    ! Loop through list, forming new one
    !----------Local vars----------!
    character(len=len(inList)+1) :: elem
    integer :: i
    logical :: myCountEmpty
    character(len=8) :: myOptions
    integer :: num
    !----------Executable part----------!
    myOptions = ' '
    if ( present(options) ) myOptions = options
    myCountEmpty = index( myOptions, 'e' ) > 0  ! .true.

    outList = inList
    if ( len_trim(inList) < 1 .or. nElement < 1 ) return
    num = NumStringElements( inList, myCountEmpty, inSeparator=inSeparator )
    if ( nElement < 1 .or. nElement > num ) return
    outList = ' '
    do i=1, num
      if ( i == nElement ) cycle
      call GetStringElement(inList, elem, i, &
        & myCountEmpty, inSeparator=inSeparator )
      outList = CatLists( outList, trim(elem), inseparator )
    enddo
  end subroutine RemoveNumFromList

  ! --------------------------------------------------  RemoveOption  -----
  subroutine RemoveOption( inOptions, outOptions, option, &
    & pattern, delims )
    ! Args
    character(len=*), intent(in)              :: inOptions
    character(len=*), intent(in)              :: option
    character(len=*), intent(out)             :: outOptions
    integer, optional, intent(in)             :: pattern
    character (len=1), dimension(2), optional, &
     & intent(in)                             :: delims
    ! Internal variables
    integer :: bloc
    logical, parameter :: COUNTEMPTY = .true.
    character (len=len(inOptions))             :: listBloc ! space-separated
    integer :: j
    integer :: myPattern
    integer :: numDashes
    ! Executable
    outOptions = inOptions
    if ( len_trim(option) < 1 ) return
    outOptions = ' '
    myPattern = 0
    if ( present(pattern) ) then
      if ( any(Pattern == (/0, 1, 2, 3, 4 /)) &  ! accept legal values only
        & ) myPattern = pattern
    endif
    if ( option(1:2) == '--' ) then
      numDashes = 2
    elseif ( option(1:1) == '-' ) then
      numDashes = 1
    else
      numDashes = 0
    endif
    select case (myPattern)
    case ( 0,1 )
      ! '-ab[arg] --xyz=arg
      ! '-a -b arg --xyz=arg'
      do bloc = 1, NumStringElements( inOptions, countEmpty, inseparator=' ' )
        listBloc = StringElement( inOptions, bloc, countEmpty, inseparator=' ' )
        ! does this block begin with one "-" or two?
        if ( numDashes == 2 .and. listBloc(1:2) == '--' ) then
          if ( listBloc == option ) listBloc = ' '
        elseif ( listBloc(1:2) == '--' ) then
          ! no match
        elseif ( numDashes == 1 .and. listBloc(1:1) == '-' ) then
          j = index( listBloc, option(2:2) )
          listBloc = listBloc(1:j-1) // listBloc(j+1:)
        endif
        outOptions = trim(outOptions) // ' ' // listBloc
      enddo
    case ( 2, 3 )
      ! '-a -b arg -xyz arg' or
      ! '-a -b arg --xyz arg'
      print *, 'Sorry--unable to remove option for pattern ', myPattern
      stop
    case ( 4 )
      ! '-a -barg --xyz arg'
      print *, 'Sorry--unable to remove option for pattern ', myPattern
      stop
    case default
      print *, 'Sorry--unable to remove option for pattern ', myPattern
      stop
    end select
  end subroutine RemoveOption

  ! --------------------------------------------------  RemoveSwitchFromList  -----
  subroutine RemoveSwitchFromList( inSwitches, outSwitches, switch, &
    & inseparator, options )
    ! Args
    character(len=*), intent(in)              :: inSwitches
    character(len=*), intent(in)              :: switch
    character(len=*), intent(out)             :: outSwitches
    character (len=*), optional, intent(in)   :: inseparator
    character (len=*), optional, intent(in)   :: options
    ! Internal variables
    integer :: i, n, details
    logical :: countEmpty
    character(len=len(inSwitches)) :: aSwitch, bareSwitch
    character(len=8) :: myOptions
    ! Executable
    myOptions = ' '
    if ( present(options) ) myOptions = options
    CountEmpty = index( myOptions, 'e' ) > 0  ! .true.
    outSwitches = ' '
    n = NumStringElements( inSwitches, countEmpty )
    if ( n < 1 ) return
    do i=1, n
      call GetStringElement( trim(inSwitches), aSwitch, i, countEmpty )
      call SplitDetails( aSwitch, bareSwitch, details )
      if ( .not. streq( Switch, bareSwitch, '-f' ) ) &
        & outSwitches = CatLists( outSwitches, aSwitch )
    enddo
  end subroutine RemoveSwitchFromList

  ! --------------------------------------------------  ReplaceSubString  -----
  subroutine ReplaceSubString (str, outstr, sub1, sub2, which, no_trim)
    ! Takes a string and replaces occurrence(s) of sub1 with sub2
    ! Defaults to replacing only the first
    ! But if which == 'all' replaces all
    ! or if which == 'last' replaces last
    ! Note that, depending on no_trim, 'all' does the following:
    ! (a) if no_trim == TRUE, multiple passes (up to 100) until no
    !     further replacements are possible
    !     ( which could be bad; e.g., if sub1 is 'sub1' and sub2 is 'sub11'
    !     then (blah)sub1(blah)sub1..' becomes '(blah)sub1111...' )
    ! (b) if no_trim is FALSE or missing, a single pass after chopping
    !     the string up into separate sub1-containing pieces
    !     ( e.g., '(blah)sub1(blah)sub1(blah)..' becomes
    !      '(blah)sub2(blah)sub2(blah)..' )
    !  
    ! Will this still work if sub1 has leading or trailing blanks? 
    ! How about sub2?
    ! do we need an optional arg, no_trim, say, that will leave them?
    ! Tried coding it, but can't say for sure it works
    !--------Argument--------!
    character (len=*), intent(in) :: str
    character (len=*), intent(in) :: sub1
    character (len=*), intent(in) :: sub2
    character (len=*) :: outstr
    character (len=*), intent(in), optional :: which
    logical, intent(in), optional :: no_trim

    !----------Local vars----------!
    integer, parameter         :: MAXREPLACEMENTS = 100
    integer :: i, array_size
    character (len=5) :: my_which
    character(len=max(len(str), len(outstr))) :: head
    character(len=max(len(str), len(outstr))) :: tail
    character(len=max(len(str), len(outstr))) :: sub_str
    character(len=max(len(str), len(outstr))), dimension(MAXREPLACEMENTS) &
      &                                      :: str_array
    logical :: my_no_trim
    !----------Executable part----------!
    head = ''
    tail = ''
    sub_str = ''
    outstr = str
    if (len_trim(str) < 1 .or. len_trim(sub1) < 1) RETURN
    my_which = 'first'
    if ( present(which) ) my_which = lowercase(which)
    my_no_trim = .false.
    if ( present(no_trim) ) my_no_trim = no_trim

    select case (my_no_trim)
    case (.false.)
      if ( index(str, trim(sub1)) < 1 ) return
      select case (my_which)
      case ('first')
        call Replace_me ( str, outstr, .false. )
      case ('last')
        call Replace_me ( str, outstr, .true. )
      case ('all')
        outstr = ' '
        str_array = ' '
        call Split_me
        do i=1, array_size
          ! print *, i, ' ', trim(str_array(i))
          call Replace_me ( trim(str_array(i)), sub_str, .false. )
          outstr = adjustl( trim(outstr) // sub_str )
        enddo
      end select
    case (.true.)
      if ( index(str, sub1) < 1 ) return
      select case (my_which)
      case ('first')
        call Replace_me_no_trim ( str, outstr, .false. )
      case ('last')
        call Replace_me_no_trim ( str, outstr, .true. )
      case ('all')
        ! Originally, I despaired of solving this
        ! CALL MLSMessage(MLSMSG_Error, ModuleName, &
        ! & 'Unable to ReplaceSubStrings with which=all and no_trim=TRUE yet')
        ! then I had an idea: Why not reinterpret this as multiple passes?
        i = 0
        str_array(1) = str
        do
          i = i + 1
          ! print *, 'i ', i
          ! print *, 'str_array(1) ', str_array(1)
          call Replace_me_no_trim ( str_array(1), str_array(2), .false. )
          ! print *, 'str_array(2) ', str_array(2)
          if ( str_array(2) == str_array(1) .or. i > MAXREPLACEMENTS ) exit
          str_array(1) = str_array(2)
        enddo
        outstr = str_array(2)
      end select
    end select
    
    contains
      subroutine Replace_me ( the_orig, after_sub, back )
        ! This replaces an instance of sub1 with sub2 in
        ! the string the_orig
        ! Either the first instance (if back == FALSE) or the last
        ! Arguments
        character(len=*), intent(in)  :: the_orig
        character(len=*), intent(inout) :: after_sub
        logical, intent(in)           :: back
        ! Local variables
        integer :: istrt1, istrt2, ihead
        if ( index(the_orig, trim(sub1)) == 0 ) then
          after_sub = the_orig
          return
        endif
        istrt1 = index(the_orig, trim(sub1), back=back)
        istrt2 = istrt1 + len_trim(sub1)
        ihead = 1
        head = ' '
        tail = ' '
        if ( istrt1 > 1 ) then
          head = the_orig(1:istrt1-1)
          ihead = istrt1 - 1
        endif
        if ( istrt2 < len_trim(the_orig)+1 ) then
          tail = the_orig(istrt2:)
        endif
        ! print *, 'len(after_sub): ', len(after_sub)
        ! print *, 'len(head): ', len(head)
        ! print *, 'len(tail): ', len(tail)
        ! print *, 'head: ', trim(head), '  ihead: ', ihead
        ! print *, 'tail: ', trim(tail), len_trim(tail)
        ! print *, 'sub1: ', trim(sub1), len_trim(sub1)
        ! print *, 'sub2: ', trim(sub2), len_trim(sub2)
        if ( sub2 /= ' ' ) then
          after_sub = adjustl(head(1:ihead) // trim(sub2) // trim(tail))
        else
          after_sub = adjustl(head(1:ihead) // trim(tail))
        endif
      end subroutine Replace_me
      subroutine Replace_me_no_trim ( the_orig, after_sub, back )
        ! This replaces an instance of sub1 with sub2 in
        ! the string the_orig -- w/o trimming leading or trailing blanks
        ! Either the first instance (if back == FALSE) or the last
        ! Arguments
        character(len=*), intent(in)  :: the_orig
        character(len=*), intent(inout) :: after_sub
        logical, intent(in)           :: back
        ! Local variables
        integer :: istrt1, istrt2, ihead
        if ( index(the_orig, sub1) == 0 ) then
          after_sub = the_orig
          return
        endif
        istrt1 = index(the_orig, sub1, back=back)
        istrt2 = istrt1 + len(sub1)
        ihead = 0
        head = ' '
        tail = ' '
        if ( istrt1 > 1 ) then
          head = the_orig(1:istrt1-1)
          ihead = istrt1 - 1
        endif
        if ( istrt2 < len(the_orig)+1 ) then
          tail = the_orig(istrt2:)
        endif
        ! Now all the possibilities:
        ! (1) the_orig = sub1
        if ( ihead == 0 .and. istrt2 > len(the_orig) ) then
          after_sub = sub2
        ! (2) the_orig = (head)sub1
        elseif ( istrt2 > len(the_orig) ) then
          after_sub = head(1:ihead) // sub2
        ! (3) the_orig = sub1(tail)
        elseif ( ihead == 0 ) then
          after_sub = sub2 // tail
        ! (4) the_orig = (head)sub1(tail)
        else
          after_sub = head(1:ihead) // sub2 // tail
        endif
      end subroutine Replace_me_no_trim
      subroutine Split_me
        ! Will this still work if some of the str_arrays end in one or more ' '?
        ! Arguments (none)
        ! Local variables
        integer :: istrt1, istrt2
        array_size = 0
        istrt2 = 0
        do
          if ( istrt2 > len_trim(str) - 1 ) return
          if ( index(str(istrt2+1:), trim(sub1)) < 1 ) then
            array_size = min(array_size+1, MAXREPLACEMENTS)
            str_array(array_size) = str(istrt2+1:)
            return
          endif
          istrt1 = istrt2 + index(str(istrt2+1:), trim(sub1))
          array_size = min(array_size+1, MAXREPLACEMENTS)
          str_array(array_size) = str(istrt2+1:istrt1 + len_trim(sub1) - 1)
          istrt2 = istrt1 + len_trim(sub1) - 1
        enddo
      end subroutine Split_me
  end subroutine ReplaceSubString

  ! --------------------------------------------------  ReverseList  -----
  function ReverseList (str, inseparator) RESULT (outstr)
    ! takes a string list, usually comma-separated,
    ! and returns one with elements in reversed order

    ! E.g., given "alpha, beta, gamma" => "gamma, beta, alpha"

    ! Limitation:
    ! No element may be longer than MAXWORDLENGTH
    !--------Argument--------!
    character (len=*), intent(in) :: str
    character (len=len(str)) :: outstr
    character (len=*), optional, intent(in)       :: inseparator

    !----------Local vars----------!
    integer :: i, istr, irev, elem, iBuf
    integer, parameter :: MAXWORDLENGTH=80
    character (len=1)               :: separator
    character (len=1), dimension(:), allocatable :: charBuf
    character (len=MAXWORDLENGTH) :: word
! Treat consecutive separators as if enclosing an empty element
    logical, parameter :: countEmpty = .TRUE.    

    !----------Executable part----------!
    if(present(inseparator)) then
      separator = inseparator
    else
      separator = COMMA
    endif

!  Special case--only one element of str
    outstr = str
    if(len(str) == 1 .OR. INDEX(str, separator) == 0) RETURN
 
! General case
    ALLOCATE(charBuf(len(str)+1), STAT=istr)
    if (istr /= 0) then
      CALL PrintMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"charBuf")
      RETURN
    endif

    outstr = ' '

! Loop over elements
    elem = 1
    iBuf=0
    do
      CALL GetStringElement(str, word, elem, countEmpty, separator)
        if(word == separator) then
          EXIT
        elseif(iBuf > len(str)) then
          EXIT
        else
          istr = MAX(len_TRIM(word), 1)
          word = Reverse(word(:istr))
        do i=1, istr
          iBuf=iBuf+1
          charBuf(iBuf) = word(i:i)
        enddo
        iBuf=iBuf+1
        charBuf(iBuf) = separator
        elem = elem+1
      endif
    enddo

    if(charBuf(iBuf) == separator) then
      iBuf = iBuf-1
    endif

    do i=1, iBuf
      irev = iBuf - i + 1
      outstr(irev:irev) = charBuf(i)
    enddo

    deallocate (charBuf)

  end function ReverseList

  ! --------------------------------------------------  ReverseStrings  -----
  subroutine ReverseStrings (str, reverse)
    ! takes an arrays of strings
    ! and returns one with elements in reversed order

    ! E.g., given (/"alpha", "beta", "gamma"/) => (/"gamma", "beta", "alpha"/)

    !--------Argument--------!
    character (len=*), dimension(:), intent(in)  :: str
    character (len=*), dimension(:), intent(out) :: reverse

    !----------Local vars----------!
    integer :: i, n
    ! Executable
    reverse = str
    n = min( size(str), size(reverse) )
    if ( n < 2 ) return
    do i=1, n
      reverse( n - i + 1 ) = str(i)
    enddo
  end subroutine ReverseStrings

  ! --------------------------------------------------  SnipList  -----
  function SnipList (str, elem, inseparator) RESULT (outstr)
    ! takes a string list, usually comma-separated,
    ! and returns one with elem removed

    ! E.g., given 
    ! str = "alpha, beta, gamma" and elem = 'gamma' => "alpha, beta"

    ! Limitation:
    ! No element may be longer than MAXWORDLENGTH
    !--------Argument--------!
    character (len=*), intent(in) :: str
    character (len=*), intent(in) :: elem
    character (len=len(str)) :: outstr
    character (len=*), optional, intent(in)       :: inseparator

    !----------Local vars----------!
    !----------Executable part----------!
!  Special case--elem not in str
    outstr = str
    if ( len_trim(elem) < 1 ) return
    if ( index(str, trim(elem) ) < 1 )  return
 
! General case
    call RemoveElemFromList( str, outstr, elem, inseparator )

  end function SnipList

  ! ---------------------------------------------  SortArray  -----

  ! This subroutine takes an array of strings
  ! and returns the array of ordered integers
  ! sorting the array; i.e., if ss[n] is the sub-string which is
  ! the n'th element, and ia[k] is the k'th element of the integer array
  ! then 
  !         {psl[ia[k]]=ss[k], k=1..n} 
  ! yields the properly sorted array
  ! (unless leftRight equals one of {"r", "R"} 
  ! in which case 
  !         {psl[k]=ss[ia[k]], k=1..n}
  ! does the job)
  ! Identical use of ia is how you would normally 
  ! sort any other arrays associated with ss
  
  ! The sorting is ordered by ascii collating sequence:
  ! "0" < "9" < "A" < "Z" < "a" < "z"
  ! unless caseSensitive is FALSE, when "0" < "9" < "A" < "a" < "Z" < "z"

  ! As an optional arg the properly sorted array is returned, too.
  ! You may safely supply the same arg for both inStrArray and sortedArray
  
  ! The optional arg options may be used to set
  ! options contains           meaning
  ! ----------------           -------
  !        c                   case insensitive
  !        s                   shorter first
  !        r                   reverse the sorting order (of both returned arrays)
  !        S                   sort as if switches
  !        L                   LeftRight is "L" (default)
  !        R                   LeftRight is "R"
  ! if the shorterFirst is TRUE, the sorting is modified
  ! so that shorter strings come first
  ! e.g., (/'abc', 'st', 'Z', '1'/) -> (/'1', 'Z', 'st', 'abc'/)
  
  ! if shorterFirst, leading spaces are always ignored
  ! otherwise they are always significant
  ! (See SortList for contrasting treatment options)
  !  (if you want them ignored, it's easy enough: create a tempArray
  !     tempArray(1:N) = adjustl(strArray(1:N))
  !   and pass it in instead)
  
  ! "Sort as if switches" is an explanation staggering in its failure
  ! to explain. What we do is to Replace each ' ' with achar(127) which
  ! has the effect of moving shorter strings from the head of the line 
  ! to the back of the line; 
  ! e.g., 'switch4' would be sorted ahead of 'sort' instead of behind it

  ! Method:
  ! The older method was removed--instead we rely on sortp which now
  ! can sort characters, too.
  
  ! Is the distinction between 'L' and 'R' sufficiently intuitive? 
  ! 'L' means that the ia[k] appears on the Left-hand side of '='
  ! 'R' means that the ia[k] appears on the Right
  subroutine SortArray( inStrArray, outIntArray, &
    & sortedArray, options )
    ! Dummy arguments
    character (len=*), dimension(:), intent(in)   :: instrarray
    integer, dimension(:), intent(out)            :: outintarray
    character (len=*), optional, intent(in)       :: options
    character (len=*), dimension(:), optional, intent(out)  &
     &                                            :: sortedarray

    ! Local variables
    logical                                :: casesensitive
    logical, parameter                     :: DeeBUG = .false.
    integer                                :: elem, nElems
    character(len=1)                       :: LeftRight
    logical                                :: reverse
    logical                                :: shorterfirst
    logical                                :: switchable
    integer, dimension(:), allocatable     :: invBinNumber 
    integer                                :: maxStrPos
    character (len=16)                     :: myOptions  
    ! integer, dimension(size(outintarray))  :: originalintarray
    integer                                :: status
    character (len=len(inStrArray)), dimension(:), allocatable    &
      &                                    :: stringArray
    character (len=len(inStrArray))        :: theString  

    ! Executable code
    myOptions = ' '
    if ( present(options) ) myOptions = options
    caseSensitive = index(myOptions, 'c' ) == 0
    reverse = index(myOptions, 'r' ) > 0
    shorterFirst = index(myOptions, 's' ) > 0
    switchable = index(myOptions, 'S' ) > 0
    leftRight = 'L'
    if ( index(myOptions, 'R' ) > 0 ) leftRight = 'R'

    nElems = size(inStrArray)
    if ( size(outIntArray) <= 0 .or. nElems <= 0 ) then
      return
    endif
    allocate (stringArray(nElems), &
     & invBinNumber(nElems), &
     & STAT=status)
    if (status /= 0) CALL PrintMessage(MLSMSG_Error, ModuleName, &
         & MLSMSG_Allocate//"stringArray, etc. in SortArray")
    outIntArray = 0
    maxStrPos = 1                ! This will hold max string LENGTH needed
    do elem = 1, nElems    
      outIntArray(elem) = 1
      if ( ShorterFirst ) then
        maxStrPos = max(maxStrPos, len_trim(adjustl(inStrArray(elem))))
      else
        maxStrPos = max(maxStrPos, len_trim(inStrArray(elem)))
      endif
    enddo                  
    if ( DEEBUG ) then
      do elem = 1, nElems    
        print *, 'Array element ', elem, ' ', trim(inStrArray(elem))
      enddo                  
    endif
    do elem = 1, nElems
      if ( Switchable ) then
        stringArray(elem) = Replace( inStrArray(elem), ' ', achar(127) )
      else
        stringArray(elem) = inStrArray(elem)
      endif
      if ( shorterFirst ) then
        ! Trickery alert!
        ! This causes shorter strings to have more leading spaces
        ! and therefore come up first when sorted
        ! (which is why we always ignore leading spaces in inStrArray)
        theString = adjustl(stringArray(elem))
        stringArray(elem) = adjustr(theString(1:maxStrPos))
      endif
    enddo
    ! Are we case-sensitive?
    if ( .not. caseSensitive ) then
      stringArray = lowercase( stringArray )
    endif

    ! Now we let asortp do the work
    call sortp( stringArray, 1, nElems, outIntArray ) ! This is OriginalIntArray
    ! Now the outIntArray is inversely related to the OriginalIntArray:
    ! (1) stringArray(OriginalIntArray(i)) = inStrArray(i)
    ! (2) stringArray(i) = inStrArray(outIntArray(i))
    ! So: inStrArray(outIntArray(OriginalIntArray(i))) =
    ! stringArray(OriginalIntArray(i)) = inStrArray(i)
    ! Or outIntArray(OriginalIntArray(i)) = i
    ! do elem = 1, nElems  
    !   outIntArray(OriginalIntArray(elem)) = elem
    ! enddo

    ! Were we asked to reverse the sorting order?
    if ( reverse ) then
      invBinNumber = outIntArray
      outIntArray = invBinNumber( NElems:1:-1 )
    endif
    
    ! Were we asked to return the sorted array, too?
    if ( present(sortedArray) ) sortedArray = inStrArray( outIntArray )
    
    ! What about left-right inversion?
    if ( LeftRight == 'L' ) then
      ! Need to 'invert' outIntArray
      invBinNumber = outIntArray
      do elem=1, nElems
        outIntArray(invBinNumber(elem)) = elem
      enddo
    endif
    deallocate (stringArray, invBinNumber, &
     & STAT=status)
    if (status /= 0) CALL PrintMessage(MLSMSG_Error, ModuleName, &
         & MLSMSG_DeAllocate//"stringArray, etc. in SortArray")

  end subroutine SortArray

  ! ---------------------------------------------  SortList  -----

  ! This subroutine takes a (usually) comma-separated string list, interprets it
  ! as a list of individual elements and returns the array of ordered integers
  ! sorting the list; i.e., if ss[n] is the sub-string which is
  ! the n'th element, and ia[k] is the k'th element of the integer array
  ! then {psl[ia[k]]=ss[k], k=1..n} yields the properly sorted list
  ! (unless the further optional arg leftRight is also supplied and equals
  ! one of {"r", "R"} in which case {psl[k]=ss[ia[k]], k=1..n})
  ! Parallel use of ia is how you would normally 
  ! sort any other arrays associated with ss
  
  ! The sorting is ordered by ascii collating sequence:
  ! "0" < "9" < "A" < "Z" < "a" < "z"
  ! unless caseSensitive is FALSE, when "0" < "9" < "A" < "a" < "Z" < "z"

  ! if countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE
  ! if TRUE, the elements would be {'a', 'b', ' ', 'd'}

  ! As an optional arg the properly sorted list is returned, too
  ! You may safely supply the same arg for both inList and sortedList
  ! As an optional arg the separator may supplied, in case it isn't comma
  ! if the optional arg ignoreLeadingSpaces is TRUE, "a, b, c" is
  ! sorted like "a,b,c"; otherwise the leading spaces make" b, c,a"

  ! Meaning of options:
  ! (see SortArray)

  ! Method:
  ! (see SortArray)
  subroutine SortList( inList, outArray, inseparator, sortedList, options )
    ! Dummy arguments
    character (len=*), intent(in)                 :: inList
    integer, dimension(:), intent(out)            :: outArray
    character (len=*), optional, intent(in)       :: inseparator
    character (len=*), optional, intent(out)      :: sortedList
    character (len=*), optional, intent(in)       :: options

    ! Local variables
    logical                                       :: countEmpty
    logical                                       :: IgnoreLeadingSpaces
    character(len=1)                              :: LeftRight
    character (len=16)                            :: myOptions  
    integer :: nElems, status, Longestlen

    character (len=1)               :: separator
    character (len=MAXSTRELEMENTLENGTH), dimension(:), allocatable    &
      &                             :: stringArray
    logical, parameter              :: DeeBUG = .false.
    ! Executable code
    if(present(inseparator)) then
      separator = inseparator
    else
      separator = COMMA
    endif
    myOptions = ' '
    if ( present(options) ) myOptions = options
    countEmpty = index(myOptions, 'e' ) /= 0
    IgnoreLeadingSpaces = index(myOptions, 'f' ) /= 0
    leftRight = 'L'
    if ( index(lowercase(myOptions), 'r' ) > 0 ) leftRight = 'R'

    if ( DEEBUG ) then
       print *, 'Entered SortList'
       print *, 'present(inseparator)?: ', present(inseparator)
       print *, 'separator: ', separator
       print *, 'string: ', trim(inList)
    endif
    if ( size(outArray) <= 0 ) return
    outArray = 0
    nElems = NumStringElements(inList, countEmpty, inseparator, Longestlen)
    if ( nElems <= 0 ) then
      return
    elseif ( Longestlen > MAXSTRELEMENTLENGTH ) then
      call PrintMessage(MLSMSG_Error, ModuleName, &
         & "Element LENGTH too long in SortList")
      return
    endif
    allocate (stringArray(nElems), STAT=status)
    if (status /= 0) CALL PrintMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"stringArray in SortList")
    call list2Array( inList, stringArray, countEmpty, inseparator, &
     & IgnoreLeadingSpaces )
    call SortArray( stringArray(1:nElems), outArray, options=options )
    if ( present(sortedList) ) then
      if ( LeftRight == 'R' ) then
        call Array2List(stringArray(1:nElems), sortedList, &
         & inseparator, outArray, leftRight='R')
      else
        call Array2List(stringArray(1:nElems), sortedList, &
         & inseparator, outArray, leftRight='L')
      endif
    endif
    deallocate (stringArray)

  end subroutine SortList

  ! ---------------------------------------------  StringElement  -----

  ! This function takes a (usually) comma-separated string list, interprets it
  ! as a list of individual elements and returns the
  ! sub-string which is the n'th element
  ! if n is too large or small, it returns blank
  ! See also GetStringElement 
  ! (which however returns separator if n too large or small)

  function StringElement(inList, nElement, countEmpty, inseparator) &
    & result(outElement)
    ! Dummy arguments
    character (len=*), intent(in)   :: inList
    integer, intent(in)         :: nElement  ! Entry number to return
    logical, intent(in)   :: countEmpty
    character (len=*), optional, intent(in)       :: inseparator
    character (len=len(inList))  :: outElement

    ! Local variables
    character (len=1)               :: separator

    ! Executable code
    if(present(inseparator)) then
      separator = inseparator
    else
      separator = COMMA
    endif
    call GetStringElement(inList, outElement, nElement, countEmpty, inseparator)
    if ( outElement == separator ) outElement = ' '
  end function StringElement

  ! ---------------------------------------------  StringElementNum  -----

  ! This function takes a (usually) comma-separated string list, interprets it
  ! as a list of individual elements, and a test string which may be an element
  ! It returns the element number of the test string in the string list
  ! or, 0 if the test string is not found
  
  ! Any leading blanks are disregarded before making the comparison;
  ! e.g., 'stare' is the same as ' stare' and is the second element of 
  ! the list 'lex, stare, decisis'
  
  ! Note: if there are multiple matches between the test string and elements
  ! of inList we return only the first
  
  ! if you want the last instead, use ReverseList on inList && subtract
  ! the answer from nElements
  
  ! This is useful because many of the hdfeos routines *inq*() return
  ! comma-separated lists
  
  ! It will be the immediate precursor function in a hash table
  ! == aka associative array == aka dictionary
  !
  ! if countEmpty is TRUE, consecutive separators, with no chars in between,
  ! are treated as enclosing an empty element
  ! Otherwise, they are treated the same as a single separator
  ! E.g., "a,b,,d" has 4 elements if countEmpty TRUE, 3 if FALSE  

  ! As an optional arg the separator may supplied, in case it isn't comma
  ! Another optional arg, part_match, returns the number of the 
  ! first element merely found in the test string; e.g.
  ! 'won, to, tree' and test 'protocol.dat' returns 2

  ! See also GetStringElement, NumStringElements

  function StringElementNum(inList, test_string, countEmpty, &
    & inseparator, part_match) result (elem)
    ! Dummy arguments
    character (len=*), intent(in)             :: inList
    character (len=*), intent(in)             :: test_string
    logical, intent(in)                       :: countEmpty
    integer                                   :: elem
    character (len=*), optional, intent(in)   :: inseparator
    logical, optional, intent(in)             :: part_match

    ! Local variables
    integer :: nElements

    character (len=MAXELEMENTLENGTH)           :: listElement
    logical ::                                    match
    ! Executable code

    nElements = NumStringElements(inList, countEmpty, inseparator)

    if(nElements <= 0) then
      elem = 0
      RETURN
    endif
    match = .false.
    if ( present(part_match) ) match = part_match

    ! Check for matches--snipping off any leading blanks
    do elem=1, nElements
      CALL GetStringElement(inList, listElement, elem, countEmpty, inseparator)
      if ( match ) then
        if (trim(listElement) /= ' ' .and. &
          & index(trim(test_string), trim(listElement)) > 0) RETURN
      else
        if(adjustl(listElement) == adjustl(test_string)) RETURN
      endif
    enddo

    elem = 0

  end function StringElementNum

  ! ---------------------------------------------  SwitchDetail  -----

  ! This function takes a (usually) comma-separated string list, interprets it
  ! as a list of individual switches, and a test switch
  ! It returns the greatest detail number of the test switch in the list
  ! or, -1 if it is not found
  
  ! As an example, say the list of switches is
  ! "abc,def,ghi2"
  ! and the test switch is "ghi"
  ! The returned value would be 2
  ! if the test switch were "abc" the returned value would be 0
  ! if the test switch were "xyz" the returned value would be -1
  ! if the test switch were "hi2" the returned value would be -1
  ! (because the start of test doesn't match the start of any list element)
  
  ! The behavior may be modified by options flag
  ! For which see comment above
  ! A special option is -R which restores an older behavior
  ! finding the Detail of the first matched string in Inlist
  ! instead of the current method which begins by sorting Inlist
  ! and then removing any duplicate switches.
  
  ! Note:
  ! By default, options automatically includes "f", for backwards compatibility
  ! If more than one switch matches the test_switch, the results are sorted
  ! and the highest match is returned.
  !   unless one of the options is "R" (see above)
  ! If the string list contains a "*" and one of the options is "w" then
  ! the test switch is automatically present
  
  function SwitchDetail( Inlist, Test_switch, Options ) result ( Detail )
    ! Dummy arguments
    character (len=*), intent(in)             :: Inlist
    character (len=*), intent(in)             :: Test_switch
    character (len=*), intent(in), optional   :: Options
    integer                                   :: Detail

    ! Local variables
    logical, parameter                        :: COUNTEMPTY = .true.
    logical                                   :: back     
    logical                                   :: dontsort 
    integer                                   :: elem     
    integer, dimension(MaxNumSwitches)        :: iarray
    integer                                   :: k   
    character (len=MAXELEMENTLENGTH)          :: listElement
    character(len=8)                          :: myOptions
    integer                                   :: nElements
    integer                                   :: startOfDetails  ! index where 
    character (len=len(test_switch)+2)        :: switch          ! the detail number 
    character (len=len(Inlist)+2)             :: Switches        ! would start
    character (len=len(Inlist)+2)             :: tempSwitches

    ! Executable code
    myOptions = '-f'
    if ( present(options) ) myOptions = trim(lowercase(options))
    detail = -1
    
    ! May return immediately under special circumstances
    ! Are either blank?
    if ( len_trim(InList) < 1 .or. len_trim(test_switch) < 1 ) return
    ! Can we tell by 1st that principles the test_switch isn't there?
    ! 1st principles means
    ! (a) no wild cards
    ! (b) comma-separated InList
    ! (c) ',switch' not found in ',switch1,switch2,..,switchn'
    if ( index(MyOptions, '*') < 1 ) then
      if ( index(MyOptions, 'c') < 1 ) then
        tempSwitches = ',' // InList
        switch = ',' // test_switch
      else
        tempSwitches = ',' // lowercase(InList)
        switch = ',' // lowercase(test_switch)
      endif
      k = index(tempSwitches, trim(switch))
      if ( k < 1 ) return
      
      ! Another short cut
      ! Available only if the switch does not appear more than once
      ! because, if it did, we would want to detect its highest Detail
      ! So, if switch appears only once, then index will always return the same
      ! k value, both when back is TRUE or FALSE
      if ( index(tempSwitches, trim(switch), back=.true.) == k ) then
        ! (d) now check if ',switch,' found in ',switch1,switch2,..,switchn,'
        tempSwitches = trim(tempSwitches) // ','
        switch = trim(switch) // ','
        if ( index(tempSwitches, trim(switch)) > 0 ) then
          ! The ',switch,' is present; 
          ! the Detail has been left unspecified, which defaults to 0
          Detail = 0
          return
        endif
      endif
    endif

    nElements = NumStringElements(inList, countEmpty)

    if ( nElements <= 0 ) Return

    back = ( index(myOptions, 'b') > 0 ) 
    if ( index(myOptions, 'c') > 0 ) then
      switch = lowercase(test_switch)
    else
      switch = test_switch
    endif
    if ( index(myOptions, 'f') > 0 ) switch = adjustl(switch)
    dontsort = ( index(myOptions, 'R') > 0 )
    if ( dontsort ) then
      Switches = InList
    else
    ! Now we want to keep only the switch with the highest details level
    ! Sort the switches to pick out the highest detail
    ! if multiple matches are found
      call sortList( CompressString(InList), iarray, ',', switches )
      tempSwitches = switches
      call GetUniqueList( tempSwitches, Switches, nElements, &
            & ignoreLeadingSpaces=.true., options='-eSL' )
    endif

   ! Check for matches
    do elem=1, nElements
      if ( back ) then
        call GetStringElement( Switches, listElement, nElements-elem+1, countEmpty )
      else
        call GetStringElement( Switches, listElement, elem, countEmpty )
      endif
      if ( index(myOptions, 'c') > 0 ) listElement = lowercase(listElement)
      if ( index(myOptions, 'f') > 0 ) listElement = adjustl(listElement)
      if ( trim(listElement) /= ' ' .and. &
          & index(trim(listElement), trim(switch)) == 1 ) then
        startOfDetails = len_trim( switch ) + 1
        if ( startOfDetails > len_trim( listElement ) ) then
          detail = 0
          exit
        endif
        ! To prevent false matches, like "walk" being matched by "walker"
        if ( isAlphabet( listElement(startOfDetails:startOfDetails) ) ) cycle
        ! Because we have sometimes allowed a "?" to be a switch
        ! (Perhaps too permissive of us)
        call ReadIntsFromChars( trim(listElement(startOfDetails:)), detail, &
          & ignore="*?")
        return
      endif
    enddo

    if ( index(myOptions, 'w') > 0 .and. index(inList, '*') > 0 ) &
     & detail = max(detail, 0)

  end function SwitchDetail

  ! ------------------------------------------------  unquote  -----
  function unquote( str, quotes, cquotes, options ) &
    & result ( outstr )
    ! function that removes a single pair of surrounding quotes from string

    ! E.g., given "Let me see." or 'Let me see.' returns
    !    Let me see.
    ! if no surrounding quotes are found, returns string unchanged; unless
    ! (1) mismatched quotes, e.g. 'Let me see." will:
    !     remove leading quote but leave trailing quote
    ! (2) a single unpaired quote found at beginning or end, will:
    !  (a) remove it if the resulting string is non-empty; or
    !  (b) return the single unpaired quote if that was the entire str
    
    ! optional arg options controls the following behaviors
    ! if options contains          meaning
    !    -----------               -------
    !         k                    strict
    !         p                    stripany
    !         r                    reverse
    !         x                    extract
    ! if strict, exceptions (1) and (2) above disregarded
    ! i.e., surrounding quotes must match, else returns string unchanged
    
    ! if stripany, any quotes, surrounding or internal,
    ! will be removed
    
    ! if extract, returns first substring surrounded by
    ! quotes; E.g., given ([a1 a2], [a3 a4]) with quotes='[' cquotes=']' returns
    !   a1 a2
    ! (This option supersedes stripany, and is automatically strict)
    
    ! if reverse, removes any quoted strings; 
    ! E.g., given 'b[a1 a2], c[a3 a4]' with quotes='[' cquotes=']' returns
    !   'b, c'
    ! (This option supersedes stripany, and is automatically strict)
    
    ! if given optional arg quotes, removes only surrounding pair:
    ! quotes[i:i] for each i=1..len[quotes]
    ! E.g., given /a\ regexp/ with quotes='/' returns
    !    a\ regexp
    
    ! if given optional args quotes & cquotes, removes only surrounding pair:
    ! quotes[i:i] on the left, cquotes[i:i] on the right, i=1..len[quotes]
    ! E.g., given [a particle] with quotes='[' cquotes=']' returns
    !    a particle
    ! (For this case, strict matching is always on)
    
    ! Useful because the parser will return quote-surrounded strings if that's
    ! how they appear in the lcf
    
    ! Calling get_string with "strip=.true." renders this unnecessary.
    ! However, you might find another use for it, especially with
    ! feature of being able to trim other, user-supplied detritus:
    ! e.g., braces, parentheses, extraneous separators
    
    ! (Aside from switches, we haven't found such a use so far;
    ! instead see more powerful ExtractSubString or ReplaceSubString)
    
    ! Note:
    ! (1) By default, if no quotes found returns input string unchanged unless
    !     extract=TRUE, in which case returns blank
    ! (2) if len(quotes) > 1, processes them in order quotes(i:i), i=1 2 ..
    !     unless extract=TRUE in which case returns after first one found
    ! (3) Perhaps extract=TRUE should be moved from here to ExtractSubString

    !--------Argument--------!
    character(len=*), intent(in) :: str
    character(len=len(str)) :: outstr
    character(len=*), intent(in), optional :: quotes
    character(len=*), intent(in), optional :: cquotes
    ! logical, intent(in), optional :: strict
    ! logical, intent(in), optional :: stripany
    ! logical, intent(in), optional :: extract
    character(len=*), intent(in), optional :: options
    !----------Local vars----------!
    character(len=len(str)) :: tmpstr
    character(len=1), parameter :: sq=''''
    character(len=1), parameter :: dq='"'
    integer :: first, last, ult, prim
    character(len=1) :: quote, cquote
    integer :: i
    logical :: mystrict
    logical :: mystripany
    logical :: myextract
    logical :: myreverse
    !----------Executable part----------!

   ult = len_trim(str)    ! Position of last non-blank char
   prim = ult - len_trim(adjustl(str)) + 1    ! Position of 1st non-blank char
   outstr=str
      
   ! LENGTH of non-blank portion of string to be trimmed must be at least 2
   if(ult-prim+1 <= 1) then
      outstr=str
      return
   endif

   myextract=.false.
   mystrict=.false.
   mystripany=.false.
   if(present(options)) then
      myextract = index(options, 'x') > 0
      mystrict = index(options, 'k') > 0
      mystripany = index(options, 'p') > 0
      myreverse = index(options, 'r') > 0
   endif
   
   mystripany = mystripany .and. (.not. myextract) .and. (.not. myreverse)
   
   ! These are initialized so that if no matching quotes found
   ! we will return    outstr = adjustl(str)
   first = prim
   last = ult

   ! trim surrounding user-supplied marks?

   if(present(quotes)) then
      if(len_trim(quotes) <= 0) then
       outstr=str
       return
      endif
      
      ! Loop over char class in string quotes
      do i=1, len_trim(quotes)
      
         quote = quotes(i:i)
         
         ! Stripany option in force?
         if ( mystripany ) then
            ! print *, 'Replacing ', quote, ' in ', trim(outstr)
            call ReplaceSubString(outstr, tmpstr, quote, '', 'all')
            outstr = tmpstr
            ! print *, trim(outstr)
            cycle
         endif

         ! Supplied with paired left and right quotes?
         if(present(cquotes)) then
            cquote=cquotes(i:i)
            mystrict=.true.
         else
            cquote=quote
         endif

        if(myreverse) then
          tmpstr = outstr
          call RemoveAnyQuotedStrings( tmpstr, outstr, quote, cquote )
          cycle
        endif
        if(myextract) then
          if ( index(str, quote) > 0 .and. index(str, cquote) > 0 ) then
            call ExtractSubString (str, outstr, quote, cquote)
            return
          endif
        elseif(mystrict) then
          if(str(prim:prim) == quote .and. str(ult:ult) == cquote) then
               outstr=str(prim+1:ult-1)
               return
          endif
      
        else
          if(str(prim:prim) == quote) then
           first=prim+1
          endif

          if(str(ult:ult) == cquote) then
             last=ult-1
          endif
        endif

      enddo
      if ( myreverse ) return

   ! Removing substring within quotes
   elseif(myreverse) then
     tmpstr = outstr
     call RemoveAnyQuotedStrings( tmpstr, outstr, sq, sq )
     tmpstr = outstr
     call RemoveAnyQuotedStrings( tmpstr, outstr, dq, dq )
     return

   ! Extracting substring within quotes
   elseif(myextract) then
     if ( index(str, sq) > 0 ) then
        call ExtractSubString (str, outstr, sq, sq)
     elseif ( index(str, dq) > 0 ) then
        call ExtractSubString (str, outstr, dq, dq)
     endif
     return

   ! insist surrounding marks match?
   elseif(mystrict) then
      if( &
      & str(prim:prim) == str(ult:ult) &
      & .and. &
        & (str(prim:prim) == sq .or. str(prim:prim) == dq) &
        & ) then
            outstr=str(prim+1:ult-1)
          else
            outstr=str
          endif
         return
      
   elseif(str(prim:prim) == sq) then
      first=prim+1
      if(str(ult:ult) == sq) then
         last=ult-1
      endif
      
   elseif(str(prim:prim) == dq) then
      first=prim+1
      if(str(ult:ult) == dq) then
         last=ult-1
      endif

   else
      first=prim
      if(str(ult:ult) == dq .or. str(ult:ult) == sq) then
         last=ult-1
      endif

   endif

   ! Still here?
   
   if ( myextract ) then
     outstr = ' '
       return
   elseif ( mystripany ) then
       return
   elseif(last >= first) then
       outstr=str(first:last)
   else
       outstr=str
   endif
      
  contains
    subroutine RemoveAnyQuotedStrings( str, out, q, cq )
      ! Args
      character(len=*), intent(in)  :: str
      character(len=*), intent(out) :: out
      character(len=1), intent(in)  :: q, cq
      ! Internal variables
      integer :: k1, k2
      character(len=len(str)) :: tmpstr
      ! Begin
      out=str
      do
        tmpstr = out
        k1 = index(tmpstr, q)
        k2 = index(tmpstr, cq)
        if ( k2 < k1+1 .or. k2 < 1 .or. k1 < 1 ) return
        if ( k1 == 1 .and. k2 == len_trim(tmpstr) ) then
          out = ' '
          return
        elseif ( k1 == 1 ) then
          out = tmpstr(k2+1:)
        elseif ( k2 == len_trim(tmpstr) ) then
          out = tmpstr(:k1-1)
        else
          out = tmpstr(:k1-1) // tmpstr(k2+1:)
        endif
      enddo
    end subroutine RemoveAnyQuotedStrings
  end function unquote

  ! ---------------------unwrap ---------------
  ! Unwrap str by replacing separators (meaning line feeds) between 2 breaks
  ! with a single break
  ! option         values                             default     
  ! inseparator --                                      ','       
  ! break       --                                      ' '       
  ! mode        -- 'soft' or 'hard'                    'hard'
  ! mode        -- don't break within quoted strings   (none)     
  function unwrap_array( strings ) result( outstr )
    ! Args
    character (len=*), dimension(:), intent(in)   :: strings
    character (len=MAXSTRLISTLENGTH)              :: outstr
    ! Internal variables
    integer                                       :: i
    ! Executable
    outstr = ' '
    if ( size(strings) < 1 ) return
    outstr = strings(1)
    do i=2, size(strings)
      outstr = trim(outstr) // ' ' // strings(i)
    enddo
  end function unwrap_array

  function unwrap_list( str, inseparator, break ) result( outstr )
    ! Args
    character (len=*), intent(in)                 :: str
    character (len=MAXSTRLISTLENGTH)              :: outstr
    character (len=*), optional, intent(in)       :: inseparator ! if not ','
    character (len=*), optional, intent(in)       :: break ! if not ' '
    ! Internal variables
    logical, parameter                            :: countEmpty = .true.
    integer                                       :: n
    character (len=256), dimension(:), pointer    :: strings => null()
    ! Executable
    outstr = ' '
    n = NumStringElements( str, countEmpty, inseparator )
    if ( n < 1 ) return
    allocate( strings(1:n) )
    call list2Array( str, strings, countEmpty, inseparator=inseparator )
    outstr = unwrap( strings )
    deallocate( strings )
  end function unwrap_list

  ! ---------------------wrap ---------------
  ! Wrap str by putting separators (meaning line feeds) between 2 breaks
  ! so no line exceeds width
  ! option         values                             default     
  ! inseparator --                                      ','       
  ! break       --                                      ' '       
  ! mode        -- 'soft' or 'hard'                    'hard'
  ! mode        -- don't break within quoted strings   (none)     
  ! dontsqueeze -- don't squeeze consecutive spaces    false (make it true?)  
  subroutine wrap_sca( str, outstr, width, &
    & inseparator, break, mode, quotes, addedLines, dontSqueeze )
    ! Args
    character (len=*), intent(in)                 :: str
    character (len=*), intent(out)                :: outstr
    integer, intent(in)                           :: width
    character (len=*), optional, intent(in)       :: inseparator ! if not ','
    character (len=*), optional, intent(in)       :: break ! if not ' '
    character (len=*), optional, intent(in)       :: mode ! if not 'hard'
    character (len=*), optional, intent(in)       :: quotes ! don't break quoted
    integer, optional, intent(out)                :: addedLines ! by wrapping
    logical, optional, intent(in)                 :: dontSqueeze
    ! Internal variables
    logical, parameter         :: countEmpty = .true.
    logical, parameter         :: deebug = .false.
    integer                    :: i
    integer                    :: j
    integer                    :: lastPos
    logical                    :: noQuotes
    integer                    :: offset
    character(len=len(outstr)) :: partstr
    character(len=1)           :: quote
    character(len=len(outstr)) :: wrpartstr
    ! Executable
    if ( present(addedLines) ) addedLines = 0
    outstr = str
    if ( len_trim(str) <= width ) return
    if ( .not. present(quotes) ) then
      call wrap_noQuotes( str, width, outstr, &
        & inseparator=inseparator, break=break, mode=mode, &
        & addedLines=addedLines, dontSqueeze=dontSqueeze )
    else
      ! We will assume that no more than one kind of quote will appear in a str
      ! if this assumption needs to be relaxed the following will be inadequate
      noQuotes = .true.
      do i=1, len_trim(quotes)
        quote = quotes(i:i)
        if ( index(str, quote) < 1 ) cycle
        ! Quoted strings require that quotes appear in pairs
        ! therefore there cannot be an odd number of them
        if ( mod( ncopies(str, quote), 2) == 1 ) cycle
        noQuotes = .false.
        lastPos = 0
        outstr = ' '
        if ( deebug ) print *, 'quote: ' // trim(quote), &
          & 'number ', NumStringElements( str, countEmpty, quote )
        ! What we'll do is to wrap only the odd-numbered elements
        do j=1, NumStringElements( str, countEmpty, quote )
          offset = lastPos
          call GetStringElement ( str, partstr, j, countEmpty, quote )
          if ( deebug ) print *, 'j ', j, ' part: ' // trim(partstr)
          if ( mod(j, 2) == 1 ) then
            call wrap_noQuotes( partstr, width, wrpartstr, &
              & inseparator=inseparator, break=break, mode=mode, offset=offset, lastPos=lastPos, &
              & addedLines=addedLines )
            if ( deebug ) print *, 'wrapped part: ' // trim(wrpartstr)
            if ( j == 1 ) then
              outstr = wrpartstr
            else
              outstr = trim(outstr) // wrpartstr
            endif
          else
            if ( deebug ) print *, '(must not wrap quoted part)'
            outstr = trim(outstr) // quote // trim(partstr) // quote
            lastPos = lastPos + len_trim(partstr) + 2
          endif
          if ( deebug ) print *, trim(outstr)
        enddo
        exit
      enddo
      ! How did we arrive here?
      if ( noQuotes ) &
        & call wrap_noQuotes( str, width, outstr, &
        & inseparator=inseparator, break=break, mode=mode, &
        & addedLines=addedLines, dontSqueeze=dontSqueeze )
    endif
  end subroutine wrap_sca

  subroutine wrap_array( str, outstrs, width, &
    & inseparator, break, mode, quotes, addedLines, dontSqueeze )
  ! Wrap str by putting each block in separate element of output array outstrs
  ! so no element exceeds width
    ! Args
    character (len=*), intent(in)                 :: str
    character (len=*), dimension(:), intent(out)  :: outstrs
    integer, intent(in)                           :: width
    character (len=*), optional, intent(in)       :: inseparator ! if not ','
    character (len=*), optional, intent(in)       :: break ! if not ' '
    character (len=*), optional, intent(in)       :: mode ! if not 'hard'
    character (len=*), optional, intent(in)       :: quotes ! don't break quoted
    integer, optional, intent(out)                :: addedLines ! by wrapping
    logical, optional, intent(in)                 :: dontSqueeze
    ! Internal variables
    ! Executable
    if ( present(addedLines) ) addedLines = 0
    outstrs(1) = str
    if ( len_trim(str) <= width ) return
    if ( .not. present(quotes) ) then
      call wrap_noQuotes( str, width, outstrs=outstrs, &
        & inseparator=inseparator, break=break, mode=mode, &
        & addedLines=addedLines, dontSqueeze=dontSqueeze )
    else
      print *, 'We cant really wrap arrays with quotes yet'
      call wrap_noQuotes( str, width, outstrs=outstrs, &
        & inseparator=inseparator, break=break, mode=mode, &
        & addedLines=addedLines, dontSqueeze=dontSqueeze )
    endif
  end subroutine wrap_array

  subroutine wrap_noQuotes( str, width, outstr, outstrs, &
    & inseparator, break, mode, offset, lastPos, addedLines, dontSqueeze )
    ! Args
    character (len=*), intent(in)                 :: str
    integer, intent(in)                           :: width
    character (len=*), optional, intent(out)                :: outstr
    character (len=*), dimension(:), optional, intent(out)  :: outstrs
    character (len=*), optional, intent(in)       :: inseparator ! if not ','
    character (len=*), optional, intent(in)       :: break ! if not ' '
    character (len=*), optional, intent(in)       :: mode ! if not 'hard', then 'soft'
    integer, optional, intent(in)                 :: offset
    integer, optional, intent(out)                :: lastPos
    integer, optional, intent(inout)              :: addedLines ! by wrapping
    logical, optional, intent(in)                 :: dontSqueeze
    ! Internal variables
    integer                                       :: dsnext    
    integer                                       :: dsp       
    integer                                       :: istr      
    integer                                       :: ko        
    integer                                       :: kp        
    character(len=1)                              :: myBreak   
    character(len=1)                              :: myMode    
    integer                                       :: myLastPos 
    integer                                       :: myOffset  
    integer                                       :: nextwidth 
    logical                                       :: NoConsecutiveSpaces
    integer                                       :: so        
    integer                                       :: sp        
    character(len=4)                              :: separator 
    ! Executable
    if(present(inseparator)) then
      separator = inseparator
    else
      separator = comma
    endif
    myBreak = ' '
    if ( present(break) ) myBreak = break
    myMode = 'h' ! hard
    if ( present(mode) ) myMode = mode
    myOffset = 0
    if ( present(offset) ) myOffset = offset
    if ( present(outstr) ) outstr = str
    if ( present(outstrs) ) outstrs(1) = str
    if ( present(addedLines) ) addedLines = 0
    NoConsecutiveSpaces = .true.
    if ( present(dontSqueeze) ) NoConsecutiveSpaces = .not. dontSqueeze
    ! print *, 'str ', trim(str)
    ! print *, 'myMode ', myMode
    ! print *, 'separator ', separator
    ! print *, 'len_trim(separator) ', len_trim(separator)
    istr = 0
    if ( len_trim(str) <= width ) return
    so = 1 ! this is the current character number of str
    ko = 1 ! this is the current character number of outstr
    myLastPos = ko + myOffset
    do
      ! how big is next width?
      nextwidth = min(width, len_trim(str) - so + 1)
      if ( so == 1 ) nextwidth = max( 1, nextwidth-myOffset )
      ! print *, 'nextwidth ', nextwidth
      if ( nextwidth < 1 ) exit
      ! does the rest of str fit within nextwidth? If so, copy it to outstr, and we're done
      if ( nextwidth >= len_trim(str) - so + 1 ) then
        if ( present(outstr) ) outstr(ko:ko + nextwidth - 1) = str(so:so + nextwidth - 1)
        if ( present(outstrs) ) then
          istr = istr + 1
          outstrs(istr) = str(so:so + nextwidth - 1)
        endif
        myLastPos = myLastPos + nextWidth
        exit
      endif
      ! Is there a separator between here and nextwidth?
      dsp = ( index( str(so:so+nextwidth-1), trim(separator), back=.true. ) )
      if ( dsp > 0 ) then
        ! Yes, so we break there
        myLastPos = 1
        sp = so + dsp - 2 + len_trim(separator)
        kp = ko + dsp - 2 + len_trim(separator)
        if ( present(outstr) ) outstr(ko:kp) = str(so:sp)
        if ( present(outstrs) ) then
          istr = istr + 1
          outstrs(istr) = str(so:sp)
        endif
        ko = kp + 1
        so = sp + 1
        if ( present(addedLines) ) addedLines = addedLines + 1
        cycle
      endif
      select case (lowerCase(myMode))
      case ('h')
        ! 'hard' wrap
        ! we will wrap to exactly width or less, even if we have to hyphenate
        ! do we have any breakable spaces in next width?
        ! dsp = index( str(so:so+nextwidth-1), ' ', back=.true. )
        dsp = index( str(so:so+nextwidth-1), myBreak, back=.true. )
        ! print *, 'so, ko ', so, ko
        ! print *, 'dsp ', dsp
        if ( dsp > 0 .and. dsp < nextwidth+1 ) then
          ! Yes, so we break there
          myLastPos = 1
          sp = so - 1 + dsp
          kp = ko + dsp - 2 + len_trim(separator)
          if ( present(outstr) ) outstr(ko:kp) = str(so:sp-1) // trim(separator)
          ko = ko + dsp + len_trim(separator) - 1
          if ( present(addedLines) ) addedLines = addedLines + 1
          ! Now treat possibility that next chars might be spaces, too
          if ( len_trim(myBreak) > 0 ) sp = sp + 1
          dsnext = findFirst( trim_safe(str(sp:)), ' ', reverse=.true. )
          if ( dsnext < 1 ) exit
          so = sp + dsnext - 1
        else
          ! No, so we hyphenate
          myLastPos = 1
          kp = ko + nextwidth - 2 + len_trim(separator)
          ! print *, 'ko, kp ', ko, kp
          ! print *, 'so, so+nextwidth-3 ', so, so+nextwidth-3
          if ( present(outstr) ) outstr(ko:kp) = str(so:so+nextwidth-3) // '-' // trim(separator)
          ko = ko + nextwidth
          so = so + nextwidth - 2
          if ( present(addedLines) ) addedLines = addedLines + 1
        endif
      case ('s')
        ! 'soft' wrap
        ! We will find the next break and wrap to that
        ! even though the resulting width may be slightly greater than planned
        ! 1st: try to wrap within width
        dsp = index( str(so:so+nextwidth-1), myBreak, back=.true. )
        ! print *, '(nextwidth+1: ', nextwidth+1
        ! print *, '(break) dsp: ', dsp
        if ( dsp > 0 .and. dsp < nextwidth+1 ) then
          myLastPos = 1
          ! Yes, so we break there
          sp = so + dsp - 1
          kp = ko + dsp - 1 + len_trim(separator)
          if ( present(outstr) ) outstr(ko:) = str(so:sp) // trim(separator)
          ! print *, 'in:  ', str(so:sp)
          ! print *, 'out: ', outstr(ko:kp)
          ko = kp + 1
          if ( present(addedLines) ) addedLines = addedLines + 1
          ! Now treat possibility that next chars might be spaces, too
          if ( len_trim(myBreak) > 0 ) sp = sp + 1
          dsnext = findFirst( trim_safe(str(sp:)), ' ', reverse=.true. )
          if ( dsnext < 1 ) exit
          sp = sp + dsnext
          so = sp - 1
        else
          ! Look for next break starting with width
          dsp = index( trim_safe(str(so+nextwidth-1:)), myBreak )
          ! print *, '(break in width) dsp: ', dsp
          if ( dsp > 0 ) then
            myLastPos = 1
            ! Yes, so we break there
            dsp = dsp + nextwidth - 1
            sp = so - 1 + dsp
            kp = ko + dsp - 2 + len_trim(separator)
            if ( present(outstr) ) outstr(ko:kp) = str(so:sp-1) // trim(separator)
            ! print *, outstr(ko:kp)
            ko = ko + dsp + len_trim(separator) - 1
            if ( present(addedLines) ) addedLines = addedLines + 1
            ! Now treat possibility that next chars might be spaces, too
            if ( len_trim(myBreak) > 0 ) sp = sp + 1
            dsnext = findFirst( trim_safe(str(sp:)), ' ', reverse=.true. )
            if ( dsnext < 1 ) exit
            sp = sp + dsnext
            so = sp - 1
          else
            ! No, so we must give up any further wrapping
            dsnext = min( len(str) - so, len(outstr) - ko )
            if ( present(outstr) ) outstr(ko:ko+dsnext) = str(so:so+dsnext)
            ! print *, outstr(ko:ko+dsnext)
            myLastPos = myLastPos + dsnext
            exit
          endif
        endif
      case default
        ! What were you thinking? 'h' or 's' are the only modes we coded
      end select
    enddo
    ! Remove consecutive spaces?
    ! print *, 'Leaving wrap_noQuotes'
    ! print *, trim(outstr)
    ! print *, 'NoConsecutiveSpaces ', NoConsecutiveSpaces
    ! print *, 'dontSqueeze ', dontSqueeze
    ! stop
    if ( len_trim(myBreak) == 0 .and. NoConsecutiveSpaces ) then
      if ( present(outstr) ) outstr = squeeze( outstr )
    endif
    if ( present(lastPos) ) lastPos = myLastPos
  end subroutine wrap_noQuotes

  ! --------------------------------------------------  WriteIntsToList  -----
  subroutine WriteIntsToList ( ints, List )
    ! Takes an array of ints and writes it as a string list
    ! E.g., given (/ 1, 2, 2, 3, 4, 5 /) returns '1 2 2 3 4 5'
    ! (Inverse of ReadIntsFromList)
    !--------Argument--------!
    character (len=*), intent(out)     :: List
    integer, dimension(:), intent(in)  :: ints
    ! Method:
    ! Use writeIntsToChars
    character(len=32), dimension(size(ints)) :: strs
    call WriteIntsToChars( ints, strs )
    call Array2List ( strs, List )
  end subroutine WriteIntsToList

!============================ Private ==============================
! ---------------------------------------  Deallocate_Index_Stack  -----
  subroutine Deallocate_Index_Stack
    ! internal variables
    integer :: stat
    ! Executable
    if ( .not. allocated(stack) ) return
    deallocate ( stack, stat=stat )
    stack_ptr = 0
  end subroutine Deallocate_Index_Stack

  subroutine prepOptions( options )
    ! Process options into separate optional args
    ! You should call this at the start of every procedure
    ! that uses options to set countEmpty, etc.
    ! Args:
    character(len=*), intent(in), optional  :: options
    ! Internal variables
    integer :: sepIndex
    character(len=16) :: myOptions
    ! Executable
    countEmpty          = .false.
    caseSensitive       = .true.
    ignoreLeadingSpaces = .false.
    separator           = ','
    myOptions = STRINGLISTOPTIONS
    if ( present(options) ) myOptions = options
    if ( len_trim(myOptions) > 0 ) then
      countEmpty          = ( index(myOptions, 'e') > 0 )
      caseSensitive       = ( index(myOptions, 'c') < 1 )
      ignoreLeadingSpaces = ( index(myOptions, 'f') > 0 )
      sepIndex = index(myOptions, 's')
      if ( sepIndex > 0 ) then
        sepIndex = sepIndex + 1
        if ( len_trim(myOptions) > sepIndex + 1 ) then
          if ( myOptions(sepIndex:sepIndex) == '{' ) &
            & separator = myOptions(sepIndex+1:sepIndex+1)
        elseif ( len_trim(myOptions) == sepIndex ) then
          separator = myOptions(sepIndex:sepIndex)
        endif
      else
        sepIndex = index(myOptions, '{')
        if ( sepIndex > 0 .and. len_trim(myOptions) > sepIndex ) &
          & separator = myOptions(sepIndex+1:sepIndex+1)
      endif
    endif
  end subroutine prepOptions

! -------------------------------------------------  Push  -----
  subroutine Push ( Index )
    ! Push the stack.  

    integer, intent(in) :: Index ! Whatever caller wants to send

    integer :: Stat
    type(Index_Stack_t), allocatable :: Temp_Stack(:)
    ! Executable

    if ( .not. allocated(stack) ) then
      ! If you allocate with lbound < 0, other stuff won't work.
      allocate ( stack(1), stat=stat )
      if ( stat /= 0 ) then
        print *, 'Unable to allocate stack'
        return
      end if
      ! call test_allocate ( stat, moduleName, 'Stack', &
      !   & ubounds=(/startingStackSize/), elementSize=storage_size(stack) / 8 )
      stack_ptr = 1
      stack(1)%index=index
    else
      ! Must increase stack size
      ! so we double it
      ! But limit number doublings to MAXDOUBLINGS
      allocate ( temp_stack(stack_ptr+1), stat=stat )
      if ( stat /= 0 ) then
        print *, 'Unable to allocate temp_stack'
        return
      end if
      temp_stack(:stack_ptr) = stack
      call move_alloc ( temp_stack, stack )
      stack_ptr = stack_ptr + 1
      stack(stack_ptr)%index =  index
    end if
  end subroutine Push

  subroutine Pop ( Frame )

    type(Index_Stack_t), intent(out) :: Frame

    double precision :: Delta ! Memory change, as accounted by Allocate_Deallocate
    logical :: HaveStack
    integer :: IDelta         ! Memory change, as accounted by system in kB
    logical :: MySilent
    logical :: MySysSize
    real :: T
    integer :: Total_used     ! Memory used in kilobytes (1024)
    character(len=10) :: Used

    ! Executable
    haveStack = allocated(stack)
    if ( haveStack ) haveStack = stack_ptr >= lbound(stack,1)

    if ( haveStack ) then
      frame = stack(stack_ptr)
      stack_ptr = stack_ptr - 1
    end if

  end subroutine Pop

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSStringLists
!=============================================================================

! $Log$
! Revision 2.91  2022/01/04 23:44:11  pwagner
! Changes to appease gfortran
!
! Revision 2.90  2021/04/29 22:52:13  pwagner
! Added WriteIntsToList
!
! Revision 2.89  2020/06/24 20:52:33  pwagner
! BooleanValue_log now respects precedence of 'and' over 'or'
!
! Revision 2.88  2020/06/09 21:55:10  pwagner
! Fix error caused by failure to Deallocate_Index_Stack
!
! Revision 2.87  2020/06/03 23:39:49  pwagner
! Improve BooleanValue_log; implemented GetMatchedParens
!
! Revision 2.86  2020/05/20 23:33:33  pwagner
! Tried to speed up SwitchDetail
!
! Revision 2.85  2019/11/11 21:17:45  pwagner
! subroutine wrap now takes optional arg dontSueeze
!
! Revision 2.84  2019/10/22 18:50:27  pwagner
! Fixed bug confusing r and R options in SortArray
!
! Revision 2.83  2019/10/21 23:18:01  pwagner
! SortArray may now reverse its sort order
!
! Revision 2.82  2019/07/09 22:59:54  pwagner
! Wrap may now put its output in an array
!
! Revision 2.81  2019/01/10 21:42:39  pwagner
! SwitchDetail returns the greatest Detail if multiple matches
!
! Revision 2.80  2018/12/11 01:21:43  pwagner
! No longer uses Printit_M
!
! Revision 2.79  2018/06/26 23:59:09  pwagner
! Dont go past end of inList in GetStringElement
!
! Revision 2.78  2017/12/12 21:22:12  pwagner
! Remove limit on character lengths in SortArray
!
! Revision 2.77  2017/12/07 22:06:12  pwagner
! Using sort_m instead of LexicalSort
!
! Revision 2.76  2017/09/25 17:24:19  pwagner
! New subroutine to RemoveOption from option string
!
! Revision 2.75  2017/08/23 16:43:48  pwagner
! Fixed bugs in SortArray; now Uses LexicalSort
!
! Revision 2.74  2017/01/25 21:12:36  pwagner
! Corrected bug in SwitchDetail; added CapitalizeList
!
! Revision 2.73  2016/12/16 21:57:09  pwagner
! Fixed a long-standing error in wrap; hopefully w/o committing new ones
!
! Revision 2.72  2016/12/14 01:23:21  pwagner
! Added unwrap
!
! Revision 2.71  2016/12/08 00:16:41  pwagner
! Added ReadNumsFromList
!
! Revision 2.70  2016/01/20 00:20:44  pwagner
! Added optional arg options to Intersection to allow wildcard matches
!
! Revision 2.69  2015/09/03 20:22:21  pwagner
! Fixed error in RemoveElemFromList
!
! Revision 2.68  2015/05/06 20:46:11  pwagner
! Repaired some error msgs
!
! Revision 2.67  2015/03/31 22:11:25  pwagner
! All args to optionDetail are optional now
!
! Revision 2.66  2014/08/19 23:15:16  vsnyder
! Added SeparatorLocation argument to GetStringElement
!
! Revision 2.65  2014/08/05 00:16:28  pwagner
! EvaluateFormula geberic: can work with Lists or Arrays
!
! Revision 2.64  2014/01/09 00:25:42  pwagner
! Added nCharsinFormat function
!
! Revision 2.63  2013/09/14 01:20:25  vsnyder
! Delete unused use name
!
! Revision 2.62  2013/09/12 23:26:47  pwagner
! Fixed bug in converting strvalues to lvalues in BooleanValue_str
!
! Revision 2.61  2013/08/28 00:38:17  pwagner
! Added a local version of PrintMessage to evade possible circular dependency
!
! Revision 2.60  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.59  2013/06/13 00:41:27  pwagner
! Removed lots of unused orts
!
! Revision 2.58  2013/05/23 16:08:20  pwagner
! Fixed bug calling unused procedure--why didn't NAG catch this?
!
! Revision 2.57  2013/05/22 20:25:44  pwagner
! Can insert, remove hash elements, scalar or array-valued
!
! Revision 2.56  2013/05/16 18:18:37  pwagner
! Corrected bugs in BooleanValue, HashElement procedures
!
! Revision 2.55  2013/05/07 21:01:21  pwagner
! Added array versions of Get, Put hash elements
!
! Revision 2.54  2013/04/05 00:47:42  pwagner
! Increased MAXSTRLISTLENGTH by factor of 4
!
! Revision 2.53  2013/04/04 22:31:05  pwagner
! Added "b"ackward option to switchDetail
!
! Revision 2.52  2012/08/30 20:51:45  pwagner
! Added RemoveSwitchFromList
!
! Revision 2.51  2012/08/27 22:54:58  pwagner
! Changed api for sorts, more useful for sorting switches
!
! Revision 2.50  2012/07/20 17:01:06  pwagner
! Added EvaluateFormula and LoopOverFormula
!
! Revision 2.49  2012/07/11 20:01:43  pwagner
! Fixed something only NAG complained about
!
! Revision 2.48  2012/07/10 15:17:15  pwagner
! Changes to GetUnique.. to work with Switches better
!
! Revision 2.47  2012/06/27 17:51:57  pwagner
! countEmpty now optional arg to remove..FromList
!
! Revision 2.46  2012/05/08 17:44:30  pwagner
! BooleanValue now a generic: values may be a string list
!
! Revision 2.45  2012/01/05 01:18:33  pwagner
! Capitalized USEd stuff
!
! Revision 2.44  2011/04/20 16:35:15  pwagner
! Added SnipList
!
! Revision 2.43  2011/02/18 18:00:10  pwagner
! Improved optionDetail; added parseOptions to fully parse commandline
!
! Revision 2.42  2010/11/05 22:23:01  pwagner
! Fixed bugs in optionDetail
!
! Revision 2.41  2010/11/03 18:29:07  pwagner
! Added optionDetail to tell whether an option is present
!
! Revision 2.40  2010/06/22 16:51:32  pwagner
! default options for switchDetail is '-f'
!
! Revision 2.39  2010/04/16 23:38:34  pwagner
! Repaired bug in switchDetail which, e.g., caused '-Sl2pc' to always return '2'
!
! Revision 2.38  2009/06/23 18:22:49  pwagner
! Added ReadIntsfromList
!
! Revision 2.37  2009/06/16 17:07:05  pwagner
! Added BuildHash to build keys, values arrays from Constructor
!
! Revision 2.36  2008/12/11 19:39:20  pwagner
! Added print statement to not_used_here
!
! Revision 2.35  2008/05/21 20:00:19  pwagner
! Must not increment optional arg unless present
!
! Revision 2.34  2008/05/09 00:24:08  pwagner
! New features added to wrap; useful for wrapLines and l2cf
!
! Revision 2.33  2008/05/02 00:08:13  pwagner
! wrap subroutine may now operate in soft mode
!
! Revision 2.32  2008/01/23 21:24:43  pwagner
! RemoveNumFromList works with inseparator correctly
!
! Revision 2.31  2007/12/19 01:28:29  pwagner
! Removed unused variables
!
! Revision 2.30  2007/09/20 17:39:59  pwagner
! Added wrap procedure
!
! Revision 2.29  2007/07/31 22:46:51  pwagner
! Added listMatches
!
! Revision 2.28  2007/06/26 00:19:21  pwagner
! Workaround another Intel bug; may expand string range into reals
!
! Revision 2.27  2007/06/21 00:49:52  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.26  2007/05/22 20:56:02  vsnyder
! don't use list-directed write to internal files
!
! Revision 2.25  2007/05/14 21:51:51  pwagner
! Bugfix for way ifc writes ints to strings
!
! Revision 2.24  2007/04/20 22:27:49  pwagner
! Minor change to keep buggy Intel compiler from producing code that bombs
!
! Revision 2.23  2006/07/12 20:37:44  pwagner
! inseparator may be any length; even 0
!
! Revision 2.22  2006/04/21 23:57:05  pwagner
! Small correction to comments on api for SwitchDetail
!
! Revision 2.21  2006/03/03 23:06:35  pwagner
! Added Intersection function
!
! Revision 2.20  2006/02/24 01:14:54  pwagner
! Added BooleanValue to evaluate boolean formulas
!
! Revision 2.19  2006/02/21 19:06:25  pwagner
! Made Get, PutHashElement routines generic
!
! Revision 2.18  2006/02/16 00:59:08  pwagner
! Fixed bug preventing "?" switch from working properly
!
! Revision 2.17  2006/01/26 00:31:46  pwagner
! Added RemoveNumFromList, MakeStringHashElement
!
! Revision 2.16  2005/11/11 21:39:12  pwagner
! added stringElement function (should we keep GetStringElement?)
!
! Revision 2.15  2005/10/18 22:52:04  pwagner
! Added IsInList function
!
! Revision 2.14  2005/09/22 23:33:58  pwagner
! date conversion procedures and functions all moved into dates module
!
! Revision 2.13  2005/09/14 22:53:26  pwagner
! Added dai_to_yyyymmdd
!
! Revision 2.12  2005/08/08 23:53:18  pwagner
! utc_format never undefined in utc_to_yyyymmdd_ints
!
! Revision 2.11  2005/08/05 16:31:07  pwagner
! Added RemoveListFromList
!
! Revision 2.10  2005/07/21 23:38:18  pwagner
! Added explanation of to-be-standard character flag options
!
! Revision 2.9  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.8  2005/06/14 18:32:25  pwagner
! Added SwitchDetail
!
! Revision 2.7  2005/03/26 00:06:54  pwagner
! Repaired RemoveElemFromList; added extract option to unquote
!
! Revision 2.6  2005/02/03 19:04:58  pwagner
! Added GetUniqueInts, utc_to_date, utc_to_time
!
! Revision 2.5  2004/10/19 22:59:08  vsnyder
! Remove USE for unused R8
!
! Revision 2.4  2004/10/13 00:51:09  vsnyder
! Move HHMMSS_value to MLSStrings
!
! Revision 2.3  2004/09/16 00:16:46  pwagner
! CatLists may cat integers onto end of stringLists
!
! Revision 2.2  2004/08/05 22:47:02  pwagner
! New interfaces to ExpandStringList for ints and logicals
!
! Revision 2.1  2004/08/04 23:17:30  pwagner
! First commit
!
