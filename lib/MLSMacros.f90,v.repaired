! Copyright 2012, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE MLSMacros              !  Macros, e..g. for the l2cf
!=============================================================================

  use HIGHOUTPUT, only: BLANKSTOCOLUMN, OUTPUTNAMEDVALUE
  use IO_STUFF, only: GET_LUN
  use MLSCOMMON, only: RANGE_T
  use MLSMESSAGEMODULE, only: MLSMSG_ERROR, &
    & MLSMSG_WARNING, MLSMESSAGE
  use MLSSTRINGLISTS, only: ARRAY2LIST, EVALUATEFORMULA, EXTRACTSUBSTRING, &
    & GETHASHELEMENT, GETSTRINGELEMENT, LIST2ARRAY, LOOPOVERFORMULA, &
    & NUMSTRINGELEMENTS, REMOVENUMFROMLIST, SWITCHDETAIL
  use OUTPUT_M, only: NEWLINE, OUTPUT
  use TOGGLES, only: SWITCHES

  implicit none
  private
  
  public :: macros, macros_t, dump_macros, expand_line, read_macros
  
! === (start of toc) ===                                                 
!     c o n t e n t s                                                    
!     - - - - - - - -                                                    

!     (data types and parameters)

! macros        The module data of stored macros
!     (subroutines and functions)
! dump_macros   Dump stored macros
! expand_line   Expand all macros in the line; 
!               the expanded result may span one or multiple lines
! read_macros   Read and store a file of macro definitions
! === (end of toc) ===

! === (start of api) ===
! dump_macros ( [int details] )
! expand_line ( char* line, char* expanded )
! read_macros ( char* filename )

! Notes:
! (1) Each macro must have a unique name composed of alphanumeric characters
!     plus underscores "-"; macro names are case-sensitive
! (2) Invoking a macro consists of prefixing its name with a "!"; the result
!     will replace the macros name, including the "!" and any arguments, 
!     with its value
! (3) The macros file is a text file with a list of macro definitions:
!     We accept 3 kinds of macros definitions:
!     (a) macro_name=value
!     (b) function(arg)="character string inluding $(arg)"
!     (c) function()="character string inluding $(n)"
!     Otherwise the macros file may contains blank lines and comments
!     Comments are preceded by the "#" character
! (4) Simple macros, without arguments are replaced their values; e.g.
!     !phaseName: phase, $
!        skipDirectWritesif=Band!BandNoOff, skipRetrievalif=Band!BandNoOff
!     becomes
!     CorePlusR4AB13: phase, $
!        skipDirectWritesif=Band13Off, skipRetrievalif=Band13Off
! (5) Formula macros, with arguments are replaced their computed values; e.g.
!     !directWriteDGG(temperature_hr,Temperature-!phase)
!     becomes
!     Label, quantity=state.temperature_hr, label='Temperature-InitPtan'      
! (6) A special loop construct using Formula macros; e.g.
!     !forall(formula,arg1,arg2,..,argn)
!     becomes
!     formula(arg1)
!       ..  ..
!     formula(argn)
! === (end of api) ===

  ! This is the type to store macros set and later expanded
  integer, parameter :: MACROSTRINGLENGTH           = 10240

  type :: macros_T
    ! Two arrays bound as a character-valued hash
    ! because the values field may contain spaces, commas, and whatnot
    ! we'll use the null as both separator and terminator
    ! This keeps track of how many macros we have so far
    integer                              :: n        = 0
    ! Each macro is defined by its key (or name)
    ! and its value
    ! In most cases the key will be some short mnemonic, containing only 
    ! alphanumeric characters and perhaps underscores
    ! e.g., "phaseName"
    ! It will be expanded by replacing an invocation, '!phaseName"
    ! by its value
    character(len=MACROSTRINGLENGTH)     :: keys     = achar(0)
    ! The value, e.g. "InitUTH"
    ! Because some values may be long lists of molecules or quantities
    ! this must be much longer than keys
    character(len=128*MACROSTRINGLENGTH) :: values   = achar(0)
  end type macros_T

  type ( macros_T ), save :: macros
  
  integer, parameter :: MAXNAMELEN  = 32
  integer, parameter :: MAXVALUELEN = 4096
  integer, parameter :: MAXN        = 128
  

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains
  ! ---------------------------------------- dump_macros ------------
  subroutine dump_macros ( details )
    ! Expand all the macros in a line
    ! Each macro to be expanded must begin with a "!" bang
    ! Ignore anything after a comment character, ";"
    ! Args
    integer, optional, intent(in)               :: details ! Don't print if < 0
    ! Internal variables
    integer                                     :: i
    integer                                     :: myDetails
    character(len=MAXNAMELEN), dimension(MAXN)  :: names
    character(len=MAXVALUELEN), dimension(MAXN) :: values
    ! Executable
    myDetails = 0
    if ( present(details) ) myDetails = details
    call outputNamedValue( 'Number of macros', macros%n )
    if ( macros%n < 1 .or. myDetails < 0 ) return
    ! Convert lists to arrays
    call outputNamedValue( 'keys', trim(macros%keys) )
    call outputNamedValue( 'values', trim(macros%values) )
    call List2Array ( macros%keys, names, &
      & countEmpty=.true., inseparator=achar(0) )
    call List2Array ( macros%values, values, &
      & countEmpty=.true., inseparator=achar(0) )
    do i=1, macros%n
      call blanksToColumn( 24 )
      call output( trim(names(i) ) )
      call blanksToColumn( 48 )
      call newLine
      call output( trim(values(i) ) )
      call newLine
    enddo
  end subroutine dump_macros

  ! ---------------------------------------- expand_line ------------
  subroutine expand_line ( line, expanded )
    ! Expand all the macros in a line
    ! Each macro to be expanded must begin with a "!" bang
    ! Ignore anything after a comment character, ";"

    ! Args
    character(len=*), intent(in)  :: LINE
    character(len=*), intent(out) :: EXPANDED
    ! Local variables
    type(Range_T) :: KRANGE   ! substring range where macro appears (incl. !)
    character(len=len(expanded)) :: arg, args, name, tline, value
    character(len=len(expanded)), dimension(MAXN) &
      &                          :: argArray, results
    integer                      :: i
    integer                      :: k
    integer                      :: macroIndex
    character(len=1)             :: macroType
    integer                      :: nargs
    ! Executable
    ! An infinite loop, exiting only when there are no more macros to expand
    expanded = line
    do
      call FindLastMacro( expanded, kRange, macroIndex, args, nargs, macroType )
      if ( kRange%Bottom < 1 ) return ! That would mean all macros have been expanded
      ! call dumpRange( kRange, 'kRange' )
      ! call output( expanded(kRange%Bottom:kRange%Top), advance='yes' )
      tline = expanded
      expanded(kRange%Bottom:) = ' '
      ! Now try to evaluate the macro
      ! call outputNamedValue( 'macroIndex', macroIndex )
      if ( macroIndex == 0 ) then
        ! A forloop
        call shift_args( args, name, macroIndex )
        call List2Array( args, argArray, countEmpty=.true. )
        ! call outputnamedValue( 'macroIndex name (after shift)', macroIndex )
        ! call outputnamedValue( 'args (after shift)', args )
        if ( macroIndex > 0 ) then
          call evaluate_macro ( macroIndex, value, &
            & argArray, nargs-1, macroType )
        else
          i = index(name, '(')
          k = index(name, ')')
          arg = name(i+1:k-1)
          call LoopOverFormula ( name, arg, &
            & argArray(1:nArgs-1), results(1:nArgs-1) )
          value = results(1)
          do i=2, nArgs
            value = trim(value) // achar(10) // results(i)
          enddo
        endif
      else
        call List2Array( args, argArray, countEmpty=.true. )
        call evaluate_macro ( macroIndex, value, argArray, nargs, macroType )
      endif
      ! call outputNamedValue( 'value', value )
      expanded(kRange%Bottom:) = trim(value) // tline(kRange%Top+1:)
      ! call output( expanded, advance='yes' )
    enddo
  end subroutine expand_line

  ! ---------------------------------------- read_macros ------------
  subroutine read_macros ( FILENAME, STATUS, MESSAGE )
    ! This reads the file of macros, filling the macros datatype

    ! Args
    character(len=*), intent(in)            :: FILENAME
    integer, intent(out), optional          :: STATUS   ! Don't use MLSMessage
    character(len=*), intent(out), optional :: MESSAGE
    ! Local variables
    integer :: I                        ! Loop inductor
    integer :: K                        ! substring position
    character(len=MACROSTRINGLENGTH) :: LINE ! A line from the file
    integer :: LUN                      ! Logical unit number
    character(len=64) :: NAME 
    integer :: NOLINES                  ! Number of lines in file
    integer :: NumMacros                ! Array size
    integer :: STAT                     ! Status flag from read
    character(len=MACROSTRINGLENGTH) :: VALUE

    ! Executable code

    ! Find a free logical unit number
    ! lun = get_lun ()
    call get_lun( lun, msg=.false. )
    if ( lun < 0 ) then
      if ( present(status) ) then
        status = -1
        if ( present(message) ) message = "No logical unit numbers available"
        return
      else
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & "No logical unit numbers available" )
      endif
    endif
    open ( unit=lun, file=filename,&
      & status='old', form='formatted', &
      & access='sequential', iostat=stat )
    if ( stat /= 0 ) then
      if ( present(status) ) then
        status = -1
        if ( present(message) ) message = "Unable to open slave file " // Filename
        return
      else
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & "Unable to open slave file " // Filename )
      endif
    endif

    ! Now read the file and count the lines
    NumMacros = 0
    noLines = 0
    firstTimeRound: do
      read ( unit=lun, fmt=*, iostat=stat ) line
      if ( stat < 0 ) exit firstTimeRound
      noLines = noLines + 1
      line = adjustl ( line )
      if ( line(1:1) /= '#' ) NumMacros = NumMacros + 1
    end do firstTimeRound

    ! Now rewind the file and read the names
    rewind ( lun )
    do i = 1, noLines
      read ( unit=lun, fmt='(a)' ) line
      line = adjustl ( line )
      if ( line(1:1) == '#' .or. len_trim(line) < 1 ) cycle
      if ( index(line, '=') < 1 ) then
        call MLSMessage ( MLSMSG_Warning, moduleName, &
          & "Could not parse macros line" // trim(line) )
        cycle
      elseif ( index(line, '()=') > 0 ) then
        ! Function definition type (c)
        k = index(line, '()=')
        name = line(:k+1)
        value = line(k+3:)
        ! call output( 'Function definition type (c): ' // trim(name) // ':' // &
        !   & trim(value), advance='yes' )
        call put_macro( name, value )
      elseif ( index(line, ')=') > 0 ) then
        ! Function definition type (b)
        k = index(line, ')=')
        name = line(:k)
        value = line(k+2:)
        ! call output( 'Function definition type (b): ' // trim(name) // ':' // &
        !   & trim(value), advance='yes' )
        call put_macro( name, value )
      else
        ! Function definition type (a)
        k = index(line, '=')
        name = line(:k-1)
        value = line(k+1:)
        call put_macro( name, value )
      endif
    end do

    close ( unit=lun )
    if ( switchDetail(switches,'macros') /=0 ) &
      & call dump_macros
  end subroutine read_macros

!--------------------------- Private Procedures ----------------------------
  subroutine dumpRange (range, name )
    ! Args
    type(Range_T), intent(in)    :: range
    character(len=*), intent(in) :: name
    ! Executable
    call outputnamedValue( 'name', name )
    call outputnamedValue( 'range', (/ range%Bottom, range%Top /) )
  end subroutine dumpRange

  subroutine evaluate_macro ( macroIndex, value, args, nargs, macroType )
    ! Args
    integer, intent(in)                         :: MACROINDEX 
    character(len=*), intent(out)               :: VALUE
    character(len=*), dimension(:), intent(in)  :: ARGS 
    integer, intent(in)                         :: NARGS
    character(len=1), intent(in)                :: MACROTYPE ! 'a', 'b', or 'c'
    ! Internal variables
    character(len=len(value))                        :: arg
    character(len=len(value)), dimension(size(args)) :: argVals, results
    character(len=len(value))                        :: formula
    integer                                          :: i, k
    character(len=len(value))                        :: name
    ! Executable
    if ( nArgs > 0 ) then
      ! call outputNamedValue ( 'number of args', nArgs )
      argVals = args
      if ( .false. ) then
      do i=1, nArgs
        call get_macro( args(i), argVals(i) )
        ! call outputnamedValue( args(i), argVals(i) )
      enddo
      endif
    endif
    call GetStringElement ( macros%keys, name, macroIndex, &
      & countEmpty=.true., inseparator=achar(0) )
    call GetStringElement ( macros%values, formula, macroIndex, &
      & countEmpty=.true., inseparator=achar(0) )
    ! call outputnamedValue( 'formula', formula )
    ! What kind of macro was it?
    select case ( macroType )
    case ( 'c' )
      ! Formula type 'c'
      ! call output( 'macro type (c): ' // trim(name), advance='yes' )
      value = EvaluateFormula ( formula, argVals(1:nArgs) )
    case ( 'b' )
      ! Formula type 'b'
      i = index(name, '(')
      k = index(name, ')')
      arg = name(i+1:k-1)
      ! call output( 'macro type (b): ' // trim(name) // ':' // trim(arg), &
      !   &  advance='yes' )
      call LoopOverFormula ( formula, arg, argVals(1:nArgs), results(1:nArgs) )
      value = results(1)
      do i=2, nArgs
        value = trim(value) // achar(10) // results(i)
      enddo
    case ( 'a' )
      ! type 'a'
      ! call output( 'macro type (a): ' // trim(name), advance='yes' )
      value = formula
    end select
  end subroutine evaluate_macro

  ! ---------------------------------------- FindLastMacro ------------
  subroutine FindLastMacro ( line, kRange, macroindex, args, nargs, macroType )
    ! Find and identify the last macro in a line
    ! Each macro must begin with a "!" bang
    ! Return its substring range kRange, its index, and any args it may have
    ! E.g., !formula( CO_HR, 400hPa, "/bin/rm -f /temp/*.txt")
    
    ! The reason we look for the last macro, is that if a formula
    ! takes an arg that happens to be a macro, too, we must
    ! expand that arg before evaluating the formula

    ! Args
    character(len=*), intent(in)  :: LINE
    type(Range_T), intent(out)    :: KRANGE
    integer, intent(out)          :: MACROINDEX
    character(len=*), intent(out) :: ARGS
    integer, intent(out)          :: NARGS
    character(len=1), intent(out) :: MACROTYPE ! 'a', 'b', or 'c'
    ! Local variables
    integer                     :: I     ! Loop inductor
    integer                     :: K     ! Substring position
    integer                     :: KF    ! Substring position
    character(len=MAXNAMELEN)   :: longername
    character(len=len(LINE))    :: tline
    ! Executable
    macroindex = -1
    kRange%Bottom = 0
    kRange%Top = 0
    args = ' '
    nargs = 0
    k = index( line, ';', back = .true. )
    if ( k < 1 .or. k > len(line) ) then
      tline = line
    else
      tline = line(1:k-1)
    endif
    k = index( tline, '!', back = .true. )
    if ( k < 1 .or. k > len(line) ) return ! Meaning there were no macros
    ! Now we assume that k marks the substring start where some
    ! stored macro name can be found. We search for that macro
    ! among the ones we have stored.
    ! Because we sort our stored macros by length, longest to shortest, 
    ! we can safely accept the first match
    do i=0, macros%n ! A trick! i==0 corresponds to the looping construct forall
      if ( i < 1 ) then
        longername = "forall("
      else
        call GetStringElement ( macros%keys, longername, i, &
          & countEmpty=.true., inseparator=achar(0) )
      endif
      macroType = 'a'
      ! Now if our macro is a formula, taking args, snip everything after
      ! the first '(' before attempting a match
      macroIndex = i
      kf = index( longername, '(' )
      if ( kf > 0 ) then
        macroType = 'b'
        if ( index( longername, '()' ) > 0 ) macroType = 'c'
        longername = longername(:kf-1)
      endif
      ! Prefix with the "!"
      longername = "!" // longername
      ! call outputNamedValue( 'longername', trim(longername) )
      ! Now is there a match?
      if ( index( tline(k:), trim(longername) ) > 0 ) exit ! A match!
    enddo
    if ( i > macros%n ) then
      ! Uh-oh, we didn't find a match
        ! call MLSMessage ( MLSMSG_Warning, moduleName, &
        !   & "Could not match macros in line" // trim(line) )
        call output( "Could not match macros in line" // trim(line), advance='yes' )
        macroIndex = -1
        return
    endif
    kRange%Bottom = k
    kRange%Top    = k + len_trim(longername) - 1 ! Because it includes '!'
    ! OK, we know name and k; now pick out args
    if ( kf < 1 ) return
    call ExtractSubString ( tline(k:), args, '(', ')' )
    if ( len_trim(args) < 1 ) return
    nargs = NumStringElements( args, countEmpty=.true., inseparator=',' )
    ! Top of kRange must point to ')' ending args
    kf = index( tline(k:), ')' )
    if ( kf < 1 ) return ! ??? How could it? Syntax error?
    kRange%Top    = k + kf - 1
  end subroutine FindLastMacro

  subroutine get_macro ( name, value )
    character(len=*), intent(in)  :: NAME 
    character(len=*), intent(out) :: VALUE
    ! Executable
    call GetHashElement ( macros%keys, macros%values, name, value, &
      & countEmpty=.true., inseparator=achar(0) )
  end subroutine get_macro

  subroutine put_macro ( name, value )
    ! We push another macro onto the end of our macros list
    ! and then sort the result so that they remain ordered
    ! from longest, by name, to shortest
    character(len=*), intent(in) :: NAME 
    character(len=*), intent(in) :: VALUE
    ! Local variables
    integer :: k
    ! Executable
    if ( macros%n < 1 ) then
      macros%n = 1
      macros%keys = trim(name) // achar(0)
      macros%values = trim(value) // achar(0)
      ! call outputNamedValue( 'keys', trim(macros%keys) )
      ! call outputNamedValue( 'values', trim(macros%values) )
      return
    endif
    k = index( macros%keys, achar(0), back=.true. )
    if ( k < len_trim(macros%keys) ) then
      macros%keys = trim(macros%keys) // achar(0) // trim(name) // achar(0)
    else
      macros%keys = trim(macros%keys) // trim(name) // achar(0)
    endif
    k = index( macros%values, achar(0), back=.true. )
    if ( k < len_trim(macros%values) ) then
      macros%values = trim(macros%values) // achar(0) // trim(value) // achar(0)
    else
      macros%values = trim(macros%values) // trim(value) // achar(0)
    endif
    macros%n = macros%n + 1
    call sort_macros
    ! call outputNamedValue( 'keys', trim(macros%keys) )
    ! call outputNamedValue( 'values', trim(macros%values) )
  end subroutine put_macro

  subroutine shift_args ( args, name, macroIndex )
    ! Shift args, returning the first as a separate "name"
    ! Args
    character(len=*), intent(inout) :: args
    character(len=*), intent(out)   :: name
    integer, intent(out)            :: macroIndex
    ! Internal variables
    integer                         :: k
    character(len=len(args))        :: targs
    ! Executable
    ! call output( 'Shifting args: ' // trim(args), advance='yes' )
    call GetStringElement ( args, name, 1, countEmpty=.true. )
    targs = args
    call RemoveNumFromList( targs, args, 1 )
    ! Now find which index number it is
    ! But we can't use StringElementNum because of the '()' nonsense
    ! macroIndex = StringElementNum( macros%keys, name, countEmpty=.true., &
    !  & inseparator=achar(0) )
    ! call outputNamedValue( 'name or formula', name )
    ! call outputNamedValue( 'remaining args', args )
    do macroIndex=1, macros%n
      call GetStringElement ( macros%keys, targs, macroIndex, &
        & countEmpty=.true., inseparator=achar(0) )
      ! call outputNamedValue( 'trying to match', targs )
      k = index(targs, '(' )
      if ( k > 0 ) then
        if ( name == targs(1:k-1) ) then
          ! call outputNamedValue( 'matched macro', targs )
          return
        endif
      endif
    enddo
    macroIndex = 0
  end subroutine shift_args

  subroutine sort_macros
    ! Sort macros by length, from longest to shortest
    ! We do this so that when we try to find a match, the
    ! first match we find will be the correct one; i.e.
    ! we won't find that !abcdef matches the macro abc
    ! Internal variables
    integer                                     :: i, j
    character(len=MAXNAMELEN)                   :: longername
    character(len=MAXVALUELEN)                  :: longervalue
    character(len=MAXNAMELEN), dimension(MAXN)  :: names
    character(len=MAXVALUELEN), dimension(MAXN) :: values
    ! Executable
    if ( macros %n < 2 ) return ! No sense sorting fewer than 2 items
    ! Convert lists to arrays
    call List2Array ( macros%keys, names, &
      & countEmpty=.true., inseparator=achar(0) )
    call List2Array ( macros%values, values, &
      & countEmpty=.true., inseparator=achar(0) )
    do i=1, macros%n-1
      do j=i+1, macros%n
        if ( len_trim(names(j)) <= len_trim(names(i)) ) cycle
        ! Uh-oh, a longer macro was found, must swap
        longername  = names(j)
        longervalue = values(j)
        names(j)  = names(i)
        values(j) = values(i)
        names(i)  = longername
        values(i) = longervalue
      enddo
    enddo
    ! Now convert from arrays back to lists
    call Array2List ( names(:macros%n),  macros%keys,   inseparator=achar(0) )
    call Array2List ( values(:macros%n), macros%values, inseparator=achar(0) )
  end subroutine sort_macros

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

END MODULE MLSMacros
!=============================================================================

!
! $Log$
! Revision 2.3  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.2  2012/08/02 23:54:21  pwagner
! Fixed bug; added optional args to bypass MLSMessage
!
! Revision 2.1  2012/08/02 17:29:23  pwagner
! First commit
!
