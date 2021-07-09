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
module MLSStrings_0               ! The lowest level string handling stuff
!=============================================================================
  implicit none
  private
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (subroutines and functions)
! Asciify            purify chars to be within printing range [32,126]
!                      (no binary) (see also ReplaceNonAscii, unAsciify)
! Capitalize         tr[a-z] -> [A-Z]
! DecimalCode        Convert binary arg to '<n>' where n is iachar(arg)
! IsAllAscii         Is a string composed entirely of ascii, i.e. non-binary
! IsAlphabet         Is the arg an alphabetical character?
! IsAscii            Is each array element ascii, i.e. non-binary
! LowerCase          tr[A-Z] -> [a-z]
! MNemonicCode       Convert binary arg to '<mnc>' 
!                      where mnc is a mnemonic code like NUL
! OctalCode          Convert binary arg to '<nnn>' 
!                      where nnn is the octal code of arg
! ReadIntsFromChars  Converts an [array of] strings to int[s] using Fortran read
! ReplaceNonAscii    Replaces every non-ascii char with newChar (see also Asciify)
! ShiftChar          shift the character so its value lies in the range 
!                     [char1,char2]
! Stretch            Insert spaces between words; optionally between letters
! Trim_safe          trims string down, but never to length 0
!
! Note: If the optional arg invert is present and TRUE,
! decimalCode, mnemonicCode, or OctalCode will perform the inverse conversion,
! e.g. from '<n>' to the binary achar(n)
! === (end of toc) ===

! === (start of api) ===
! char* Asciify ( char* str, [char* how] )
! char* Capitalize ( char* str )
! log isAllAscii( char* arg )
! log isAlphabet( char arg, [char inputCase] )
! log isAscii( char arg )
! char* LowerCase ( char* str )
! int nCharsinFormat ( char* Format )
! char(5) decimalCode( char* arg, [log invert] )
! char(5) mnemonicCode( char* arg, [log invert] )
! char(5) octalCode( char* arg, [log invert] )
! readIntsFromChars ( char* strs[(:)], int ints[(:)], &
!       & [char* forbiddens], [char* ignore] )
! char* ReplaceNonAscii ( char* str, char newChar, [char* exceptions] )
! char* stretch ( char* str, [char* options] )
! char* stretch ( char* str, int how_many, [char* options] )
! char* trim_safe ( char* str )
! === (end of api) ===

  public :: Asciify, Capitalize, DecimalCode, &
    & IsAllAscii, IsAlphabet, IsAscii, &
    & NCharsInFormat, Lowercase, &
    & MNemonicCode, OctalCode, &
    & Readintsfromchars, ReplaceNonAscii, ShiftChar, Stretch, Trim_Safe

  interface Asciify
    module procedure Asciify_Scalar, Asciify_1d, Asciify_2d, Asciify_3d
  end interface

  interface ReadIntsFromChars
    module procedure ReadAnIntFromChars, ReadIntArrayFromChars
  end interface

  interface Stretch
    module procedure Stretch_one, Stretch_many
  end interface

  ! strings2Ints
  integer, public, parameter      :: Lenorsizetoosmall=-999
  ! readAnIntFromChars
  integer, public, parameter      :: Stringcontainsforbiddens=-999

  ! The following array is used to encode non-ascii characters mnemonically
  character(len=*), dimension(33), private, parameter :: MnemonicCodes = (/ &
   & 'nul', &
   & 'soh', &
   & 'stx', &
   & 'etx', &
   & 'eot', &
   & 'enq', &
   & 'ack', &
   & 'bel', &
   & 'bs ', &
   & 'ht ', &
   & '1f ', &
   & 'vt ', &
   & 'ff ', &
   & 'cr ', &
   & 'so ', &
   & 'si ', &
   & 'dle', &
   & 'dcl', &
   & 'dc2', &
   & 'dc3', &
   & 'dc4', &
   & 'nak', &
   & 'syn', &
   & 'etb', &
   & 'can', &
   & 'em ', &
   & 'sub', &
   & 'esc', &
   & 'fs ', &
   & 'gs ', &
   & 'rs ', &
   & 'us ', &
   & 'del' /)
  
contains

  ! -------------------------------------------------  Asciify  -----
  ! takes input string and replaces any non-printing characters
  ! with corresponding ones in range [32,126]
  ! leaving other chars alone
  !
  ! How the replacement is done is according to the optional arg
  !    how         meaning
  !    ---         -------
  !  'shift'       shift to the character with matching modulus (default)
  !                  (a poor choice for a default, in my opinion)
  !  'snip'        remove the offending character (shortening the string)
  !  'decimal'     return <nnn> where nnn is the decimal value (e.g. 0)
  !  'octal'       return <nnn> where nnn is the octal value (e.g. 000)
  !  'mnemonic'    return <ID> where ID is the mnemonic code (e.g. NUL)
  !  '*'           replace with whatever character how was (e.g., '*')
  ! Note that some of these options may output a longer string than the input

  ! (see also ReplaceNonAscii, unAsciify)
  ! Note that only the effects of decimal, octal, and mnemonic may be inverted
  ! by unAsciify
  function Asciify_scalar (STR, HOW) result (OUTSTR)
    !--------Argument--------!
    character (len=*), intent(in)           :: STR
    character (len=5*len(str))              :: OUTSTR
    character (len=*), optional, intent(in) :: HOW

    !----------Local vars----------!
    integer :: I, K
    character(len=5) :: insert
    character(len=1), dimension(len(str)) :: mold
    character(len=8) :: myHow
    !----------Executable part----------!
    outstr=str
    myHow = 'shift'
    if ( present(how) ) myHow = how
    mold = transfer(str,mold,size=len(str))
    if ( all(isAscii(mold)) ) return
    outstr = ' '
    select case (myHow)
    case ('shift')
      do i=1, len(str)
        if ( isAscii(str(i:i)) ) cycle
        outstr(i:i) = ShiftChar(str(i:i))
      end do
    case ('snip')
      k = 1
      do i=1, len(str)
        if ( isAscii(str(i:i)) ) then
          outstr(k:k) = str(i:i)
          k = k + 1
        endif
      end do
    case default
    ! case ('decimal', 'octal', 'mnemonic')
      k = 1
      do i=1, len(str)
        if ( isAscii(str(i:i)) ) then
          outstr(k:k) = str(i:i)
          k = k + 1
        else
          if ( myHow == 'decimal' ) then
            insert = decimalCode( str(i:i) )
          elseif ( myHow == 'octal' ) then
            insert = octalCode( str(i:i) )
          elseif ( myHow == 'mnemonic' ) then
            insert = mnemonicCode( str(i:i) )
          else
            insert = myHow(1:1)
          endif
          outstr(k:) = insert
          k = k + len_trim(insert)
        endif
      end do
    end select
  end function Asciify_scalar

  function Asciify_1d (STR, HOW) result (OUTSTR)
    character (len=*), dimension(:), intent(in)           :: STR
    character (len=5*len(str)), dimension(size(str))      :: OUTSTR
    character (len=*), optional, intent(in) :: HOW
    integer :: i
    do i=1, size(str)
      outstr(i) = Asciify( str(i), how )
    enddo
  end function Asciify_1d

  function Asciify_2d (STR, HOW) result (OUTSTR)
    character (len=*), dimension(:,:), intent(in)                  :: STR
    character (len=5*len(str)), dimension(size(str,1),size(str,2)) :: OUTSTR
    character (len=*), optional, intent(in) :: HOW
    integer :: i
    do i=1, size(str,2)
      outstr(:,i) = Asciify( str(:,i), how )
    enddo
  end function Asciify_2d

  function Asciify_3d (STR, HOW) result (OUTSTR)
    character (len=*), dimension(:,:,:), intent(in)                  :: STR
    character (len=5*len(str)), dimension(size(str,1),size(str,2),size(str,3))&
      &  :: OUTSTR
    character (len=*), optional, intent(in) :: HOW
    integer :: i
    do i=1, size(str,3)
      outstr(:,:,i) = Asciify( str(:,:,i), how )
    enddo
  end function Asciify_3d

  ! -------------------------------------------------  CAPITALIZE  -----
  elemental function Capitalize (STR) result (OUTSTR)
    ! takes a-z and replaces with A-Z 
    ! leaving other chars alone
    !--------Argument--------!
    character (len=*), intent(in) :: STR
    character (len=len(str)) :: OUTSTR

    !----------Local vars----------!
    integer :: I, ICODE
    integer, parameter :: OFFSET=iachar("A")-iachar("a")
    !----------Executable part----------!
    outstr=str

    do i=1, len_trim(str) ! Won't need to Capitalize trailing spaces
       icode=iachar(outstr(i:i))
       if ( icode >=iachar("a") .and. icode <= iachar("z")) then
          outstr(i:i)=achar(icode+offset)
       end if
    end do

  end function Capitalize
  
  ! ---------------------------------------------------  isAllAscii  -----
  elemental function isAllAscii(arg) result(itIs)
    ! Returns TRUE if each substring of arg is in range of printing chars [32,126]
    ! Args
    character(len=*), intent(in) :: arg
    logical                      :: itIs
    ! Internal variables
    integer :: i
    ! Executable
    itIs = isAscii(arg(1:1))
    if ( len(arg) < 2 .or. .not. itIs ) return
    do i=2, len(arg)
      itis = itIs .and. isAscii(arg(i:i))
    enddo
  end function isAllAscii

  !
  ! ---------------------------------------------------  isAlphabet  -----
  elemental function isAlphabet(arg, inputcase) result(itIs)
    ! Returns TRUE if arg alphabetical; 
    ! i.e.is one of {'a', 'b', ..}
    ! Note: to check if input is UPPER  lower, either, set
    ! inputcase          case
    ! ----------         ----
    !   UPPER             u
    !   lower             l
    !   either            e (default)
    ! Args
    character(len=1), intent(in) :: arg
    character(len=1), optional, intent(in) :: inputcase
    logical                      :: itIs
    ! Internal variables
    logical :: itsEither
    logical :: itsLower
    character(len=*), parameter :: list='abcdefghijklmnopqrstuvwxyz'
    character(len=1)            :: myCase
    ! Executable
    myCase = 'e'
    if ( present(inputcase) ) myCase = inputcase
    itsLower = ( index(list, arg) > 0 )
    itsEither = ( index(list, lowercase(arg)) > 0 )
    select case(myCase)
    case ('u')
      itIs = itsEither .and. .not. itsLower
    case ('l')
      itIs = itsLower
    case default
      itis = itsEither
    end select    
  end function isAlphabet

  ! ---------------------------------------------------  isAscii  -----
  elemental function isAscii(arg) result(itIs)
    ! Returns TRUE if arg is in range of printing chars [32,126]
    ! Args
    character(len=1), intent(in) :: arg
    logical                      :: itIs
    ! Internal variables
    integer, parameter :: pcMin = iachar(' ')
    integer, parameter :: pcMax = iachar('~')
    integer :: icode
    ! Executable
    icode = iachar(arg)
    itis = .not. ( icode < pcMin .or. icode > pcMax )
  end function isAscii

  
  ! --------------------------------------------------  LowerCase  -----
  elemental function LowerCase (str) result (outstr)
    ! takes A-Z  and replaces with a-z
    ! leaving other chars alone
    !--------Argument--------!
    character (len=*), intent(in) :: STR
    character (len=len(str))      :: OUTSTR

    !----------Local vars----------!
    integer            :: i, icode
    integer, parameter :: offset=IACHAR("a")-IACHAR("A")
    !----------Executable part----------!
    outstr=str

    DO i = 1, len_trim(str) ! Won't need to Lowercase trailing spaces
       icode=IACHAR(outstr(i:i))
       IF ( icode >=IACHAR("A") .AND. icode <= IACHAR("Z")) THEN
          outstr(i:i)=achar(icode+offset)
       END IF
    END DO

  end function LowerCase

  ! ---------------------------------------------------  nCharsinFormat  -----
  function nCharsinFormat ( Format ) result( nplusm )
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
    character(len=20) :: kChar                               
    integer :: n, m, p, q
    ! Executable                                                          
    nplusm = 0                                                            
    kChar=lowerCase(Format)
    if ( index( KChar, 'es' ) > 0 ) then
      p = index( KChar, 'es' )
      kChar = kChar(:p-1) // 'f' // KChar(p+2:)
    elseif ( index( KChar, 'en' ) > 0 ) then
      p = index( KChar, 'en' )
      kChar = kChar(:p-1) // 'f' // KChar(p+2:)
    elseif ( index( KChar, 'g' ) > 0 ) then
      p = index( KChar, 'g' )
      kChar = kChar(:p-1) // 'f' // KChar(p+1:)
    elseif ( index( KChar, 'e' ) > 0 ) then
      p = index( KChar, 'e' )
      kChar = kChar(:p-1) // 'f' // KChar(p+1:)
    elseif ( index( KChar, 'd' ) > 0 ) then
      p = index( KChar, 'd' )
      kChar = kChar(:p-1) // 'f' // KChar(p+1:)
    endif
    p = index( KChar, 'f' )
    if ( p < 1 ) return
    call readIntsFromChars ( KChar(p+1:), m, ignore='.,*' )
    if (m < 1) then
      print *, 'Format: ', trim(Format)
      call PrintMessage ( ModuleName, &              
      & 'Bad conversion to m in OUTPUT_xxxLE (format not "{defg}"' )      
    endif
    if ( index(TRIM(Format), 'x' ) == 0 ) then                          
      n = 0                                                               
    else                                                                  
      p = index( KChar, '(' )
      q = index( KChar, 'x' )
      call readIntsFromChars ( KChar(p+1:q-1), n, ignore='.,*' )
      if (n < 1) then                                                     
        print *, 'Format: ', trim(Format)
        print *, trim(kChar)                                              
        print *, trim(Format)                                           
        call PrintMessage ( ModuleName, &                     
          & 'Bad conversion to n in OUTPUT_xxxLE (format not "(nx)"' )  
      end if                                                              
    end if                                                                 
    nplusm = n + m                                                        
  end function nCharsinFormat

  ! ----------------  readIntsFromChars  ----- readNumsFromChars  -----
  ! This family of routines reads either a single number or an array
  ! depending on the shape of its args
  ! If called with type integer (real, double) return args, it will try to read
  ! the string data as integers (reals, doubles)
  
  ! Use:
  ! Depending on optional arguments you can numerical part from a combination
  ! of number and unit, e.g. 6.3km
  
  ! Beware of mixing e-type format with dimensions
  ! E.g., '3.2e0ppmv' will very likely confuse the subroutine
  ! (is the unit string 'ppmv' or is it 'e0ppmv'?)
  subroutine readAnIntFromChars (str, num, forbiddens, ignore)
    ! takes a string and returns an integer
    ! using Fortran "read"
    ! (which could cause an io error--that's why this
    ! subroutine exists)
    ! If the string is blank or contains one of forbiddens
    ! the int is STRINGCONTAINSFORBIDDENS
    
    ! Then snip away any from the set ignore before the first digit
    ! if ignore is present.
    ! If ignore is '*', that means ignore all alphabetical chars
    ! If ignore contains '*', that means ignore all alphabetical chars
    ! plus any other chars among ignore.
    ! Note that '*' will only escape alphabetical chars.
    ! If you want to escape any other chars in addition to alphabeticals,
    ! you should add that char to the escape string 
    ! (e.g. ':*' to escape all alphabeticals and to escape ':', too)
    ! If the string is composed entirely of ignorable chars, int is 0
    ! If the string contains multiple numbers, separated by ignorables.
    ! only the first number is returned.
    
    ! Finally attempt to read as an int what remains
    ! If that should fail as a last resort return STRINGCONTAINSFORBIDDENS

    ! Examples:
    ! (1) if str='band13a' and ignore='*', int will be 13
    ! (2) if str='3 cm' and forbiddens='c', int will be left undefined
    !     because of the 'm'
    ! (3) if str='b7f2' and ignore='*', int will be 7

    ! Limitation: you're unable to "escape" a * so you'll have to
    ! preprocess the * away if you really want to read a string which has
    ! a * in it somewhere
    !
    !--------Argument--------!
    character (len=*), intent(in) ::   str
    integer, intent(out)          ::   num
    include 'ReadANumFromChars.f9h'
  end subroutine readAnIntFromChars

  subroutine readIntArrayFromChars (strs, ints, forbiddens, ignore)
    ! takes an array of strings and returns integer array
    ! using Fortran "read"
    ! If any element of string array is blank or contains one of forbiddens
    ! the corresponding element of ints is left undefined
    ! Not useful yet
    !
    !--------Argument--------!
    !    dimensions are (len(strs(1)), size(strs(:)))
    character (len=*), intent(in), dimension(:) ::   strs
    integer, intent(out), dimension(:)          ::   ints
    character (len=*), intent(in), optional     ::   forbiddens
    character (len=*), optional, intent(in)     ::   ignore

    !----------Local vars----------!
    integer :: i, arrSize
    !----------Executable part----------!

   ! Check that all is well (if not returns blanks)
   arrSize = MIN(size(strs), size(ints))
   if ( arrSize <= 0 ) then
     ints = LENORSIZETOOSMALL
     return
   endif
   do i=1, arrSize
     call readAnIntFromChars(strs(i), ints(i), forbiddens, ignore)
   enddo

  end subroutine readIntArrayFromChars

   ! --------------------------------------------------  ReplaceNonAscii  -----
  function ReplaceNonAscii (str, newchar, exceptions) result (outstr)
    ! takes a string and returns one with non-ascii chars replaced by newChar
    ! E.g., to replace every char(0), which is the NUL character, 
    ! and a trailing char(13), which is a line feed,
    ! with blanks, which is char(32)
    ! arg = ReplaceNonAscii( arg, char(32) )
    
    ! If exceptions are present, don't replace them
    
    ! (see also Asciify)
    character(len=*), intent(in)           :: str
    character(len=1), intent(in)           :: newChar
    character(len=*), intent(in), optional :: exceptions
    character(len=len(str))                :: outstr
    ! Internal variables
    integer :: i
    ! Executable
    outstr = str
    do i=1, len(str)
      if ( present(exceptions) ) then
        if ( index(exceptions, str(i:i)) > 0 ) cycle
      endif
      if ( .not. isAscii(str(i:i)) ) outstr(i:i) = newChar
    enddo
  end function ReplaceNonAscii

  ! ------------------------------------------------- stretch --------
  ! Insert extra spaces between words
  ! E.g., turns 'How long, America' into 'How  long,  America'
  ! optionally inserts a space after every character
  ! making it 'H o w   l o n g ,   A m e r i c a'
  
  ! options, if present, can contain the following characters
  !  character                 effect
  !     a                   insert space after every character
  !    o[xyz..]             insert space only after x or y or z or ..
  function stretch_many( str, how_many, options ) result( stretched )
    ! Insert how_many spaces instead of just one
    ! Args
    character(len=*), intent(in)              :: str
    integer, intent(in)                       :: how_many
    character(len=*), optional, intent(in)    :: options
    character(len=(how_many+1)*(len(str)+1))  :: stretched
    ! Internal variables
    integer                                   :: i
    ! Executable
    stretched  = str
    do i=1, how_many
      stretched = stretch_one( stretched, options )
    enddo
  end function stretch_many

  function stretch_one( str, options ) result( stretched )
    ! Args
    character(len=*), intent(in)           :: str
    character(len=*), optional, intent(in) :: options
    character(len=2*len(str)+1)            :: stretched
    ! Internal variables
    integer                                :: cpos ! current position in str
    integer                                :: cposq ! current position in stretched
    logical                                :: newWord
    logical                                :: everywhere
    character(len=1)                       :: space
    character(len=256)                     :: onlyAfter
    ! Executable
    everywhere = .false.
    if ( present(options) ) everywhere = ( index(options, 'a') > 0 )
    onlyAfter = ' '
    if ( present(options) ) then
      cpos = index(options, 'o')
      if ( cpos > 0 ) then
        cpos = index(options, '[')
        cposq = index(options, ']')
        onlyAfter = options(cpos:cposq)
      endif
    endif
    space = ' '
    stretched = str
    if ( len_trim(str) < 2 ) return
    stretched = ' '
    if ( everywhere ) then
      cpos = len_trim(str)
    ! Fortran does not support (start : end : stride) syntax for substrings
    ! stretched(1:2*cpos+1:2) = str(1:cpos)
    ! stretched(2:2*cpos:2) = ' '
      do cpos = 1, len_trim(str)
        ! This is easy -- snip every space no matter where
        stretched(2*cpos-1:2*cpos) = str(cpos:cpos) // ' '
      enddo
    elseif( len_trim(onlyAfter) > 0 ) then
      stretched = str
      cposq = 1
      if ( len_trim(str) < 2 ) return
      do cpos = 2, len_trim(str)
        cposq = cposq + 1
        stretched(cposq:cposq) = str(cpos:cpos)
        if ( cpos == len_trim(str) ) then
          cposq = cposq + 1
          stretched(cposq:cposq) = ' '
        elseif ( index(trim(onlyAfter), str(cpos:cpos)) > 0 .and. &
          & str(cpos+1:cpos+1) /= ' ' ) then
          cposq = cposq + 1
          stretched(cposq:cposq) = ' '
        endif
      enddo
    else
      stretched(1:1) = str(1:1)
      cposq = 1
      newWord = ( str(1:1) == space )
      do cpos = 2, len_trim(str)
        if ( newWord ) then
          ! Already have at least one space, so add one more
          cposq = cposq + 1
          stretched(cposq:cposq) = ' '
          if ( str(cpos:cpos) /= space ) then
            cposq = cposq + 1
            stretched(cposq:cposq) = str(cpos:cpos)
          endif
        else
          ! don't insert
          cposq = cposq + 1
          stretched(cposq:cposq) = str(cpos:cpos)
        endif
        ! Have we reached a space which divides words?
        newWord = ( str(cpos:cpos) == space )
      enddo
    endif
  end function stretch_one
       
  ! -------------------------------------------------  TRIM_SAFE  -----
  function trim_safe (STR) result (OUTSTR)
    ! trims str returning a string of length no less than 1
    ! similar to trim, but will return a single blank character
    ! Useful in those cases where trim would result in strings of length 0
    ! E.g., MakeHDFAttribute(trim(' ')) fails but
    ! E.g., MakeHDFAttribute(trim_safe(' ')) succeeds
    !--------Argument--------!
    character (len=*), intent(in) :: STR
    character (len=max(len_trim(str), 1)) :: OUTSTR

    !----------Executable part----------!
    outstr=' '

    if ( len_trim(str) > 0 ) outstr=trim(str)

  end function trim_safe

  ! ---------------------------------------------------  octalCode  -----
  function octalCode(arg, invert) result(theCode)
    character(len=*), intent(in)  :: arg
    logical, optional, intent(in) :: invert
    character(len=5) :: theCode
    ! Returns '<nnn>' where nnn is the octal code of arg
    integer :: iCode
    character(len=3) :: nnn
    logical :: myInvert
    ! Executable
    ! print *, 'octal arg ', arg
    myInvert = .false.
    if ( present(invert) ) myInvert = invert
    if ( myInvert ) then
      read( arg(2:len_trim(arg)-1),'(o3.3)') iCode
      ! print *, 'arg(2:len_trim(arg)-1) ', arg(2:len_trim(arg)-1)
      ! print *, 'iCode ', iCode
      theCode = achar(iCode)
    else
      write(nnn,'(o3.3)') iachar(arg)
      theCode = '<' // nnn // '>'
    endif
  end function octalCode

  ! ---------------------------------------------------  decimalCode  -----
  function decimalCode(arg, invert) result(theCode)
    character(len=*), intent(in)  :: arg
    logical, optional, intent(in) :: invert
    character(len=5) :: theCode
    ! Returns '<n>' where n is iachar(arg)
    integer :: iCode
    character(len=3) :: n
    logical :: myInvert
    ! Executable
    ! print *, 'decimal arg ', arg
    myInvert = .false.
    if ( present(invert) ) myInvert = invert
    if ( myInvert ) then
      call readIntsFromChars( arg(2:len_trim(arg)-1), iCode )
      theCode = achar(iCode)
      ! print *, 'arg(2:len_trim(arg)-1) ', arg(2:len_trim(arg)-1)
      ! print *, 'iCode ', iCode
    else
      write(n,'(i3)') iachar(arg)
      theCode = '<' // trim(adjustl(n)) // '>'
    endif
  end function decimalCode

  ! ---------------------------------------------------  mnemonicCode  -----
  function mnemonicCode(arg, invert) result(theCode)
    character(len=*), intent(in) :: arg
    logical, optional, intent(in) :: invert
    character(len=5) :: theCode
    ! Returns '<mnc>' where mnc is a mnemonic code like NUL
    integer :: icode
    logical :: myInvert
    ! Executable
    ! print *, 'mnemonic ', arg
    myInvert = .false.
    if ( present(invert) ) myInvert = invert
    if ( myInvert ) then
      ! Avoiding FindFirst so this module stands alone
      do iCode=1, size(MNEMONICCODES)
        if ( MNEMONICCODES(iCode) == lowercase(arg(2:len_trim(arg)-1)) ) exit
      enddo
      ! iCode = FindFirst( MNEMONICCODES == lowercase(arg(2:len_trim(arg)-1)) )
      if ( iCode > size(MNEMONICCODES) ) iCode = 128
      theCode = achar(max(iCode - 1, 0))
      ! print *, 'arg(2:len_trim(arg)-1) ', arg(2:len_trim(arg)-1)
      ! print *, 'iCode ', iCode
    else
      theCode = arg
      if ( isAscii(arg) ) return
      icode = iachar(arg) + 1
      icode = min(icode, size(MNEMONICCODES))
      theCode = '<' // trim( MNEMONICCODES(icode) ) // '>'
    endif
  end function mnemonicCode

  ! ---------------------------------------------------  shiftChar  -----
  elemental character function shiftChar( arg, char1, char2 )
    character(len=1), intent(in)           :: arg
    character(len=1), optional, intent(in) :: char1
    character(len=1), optional, intent(in) :: char2
    ! shift the character so its value lies in the range [char1,char2]
    ! (by default char1 = ' ', char2 = '~')
    integer :: m1, m2
    integer :: icode
    integer :: resultcode
    ! Executable
    m1 = iachar(' ')
    if ( present(char1) ) m1 = iachar(char1)
    m2 = iachar('~')
    if ( present(char2) ) m2 = iachar(char2)
    icode = iachar(arg)
    shiftChar = arg
    resultcode = -1 ! This means the character did not need shifting
    if ( icode < m1 ) then
      resultcode = m1 + mod( icode+256, m2-m1 )
    elseif ( icode > m2 ) then
      resultcode = m1 + mod( icode, m2-m1 )
    endif
    if ( resultcode > -1 ) shiftChar = achar(resultcode)
  end function shiftChar

  ! ------------------------------------  PrintMessage  -----
  subroutine PrintMessage ( name, line, advance )
    ! Args
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: line
    character (len=*), intent(in), optional :: Advance ! Do not advance
    !                                 if present and the first character is 'N'
    !                                 or 'n'
    ! Local variables
    integer :: nChars
    character(len=len(line) + len(name) + 3) :: thus
    ! Executable
    nChars = len(line)
    thus = line
    if ( len_trim(name) > 0 ) then
      nChars = len(line) + len(name) + 3
      thus = '(' // trim(name) // ') ' // line
    endif
    print *, 'thus(1:nChars): ' 
    print *, thus(1:nChars)
    stop
  end subroutine PrintMessage

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSStrings_0
!=============================================================================

! $Log$
! Revision 2.3  2021/07/09 21:49:09  pwagner
! Clearer comments in Asciify
!
! Revision 2.2  2019/07/17 20:11:49  pwagner
! Stretch may now take an iteration arg, how_many
!
! Revision 2.1  2019/04/09 20:35:47  pwagner
! First commit
!
