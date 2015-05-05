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
module ManipulationUtils        ! operations to manipulate quantities
  !=============================================================================

  use Allocate_deallocate, only: allocate_test, deallocate_test
  use Dates_module, only: tai93s2hid
  use Dump_0, only: dump
  use HighOutput, only: outputNamedValue
  use MLSKinds, only: rv
  use MLSL2Options, only: MLSMessage
  use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Warning
  use MLSFinds, only: findFirst, findLast, findNext
  use MLSStats1, only: MLSCount, MLSMin, MLSMax, MLSMean, MLSMedian, &
    & MLSRMS, MLSstddev
  use MLSStringLists, only: array2List, catLists, getStringElement, &
    & List2Array, numStringElements, &
    & ReplaceSubstring
  use MLSStrings, only: indexes, lowerCase, nCopies, splitNest, stretch
  use Output_m, only: output
  use VectorsModule, only: vectorValue_t, m_fill, reshapeVectorValue
  ! This module allows us to do algebraic operations on vector quantities
  ! saving the result in a vector quantity
  ! See also Algebra Module (though I never got Algebra to work--paw)

  implicit none
  private
  public :: manipulate

! === (start of toc) ===
! Manipulate     Apply manipulation encoded in a string m to fill quantity
!                     using other quantities a, b, and possibly constant c
! === (end of toc) ===

! === (start of api) ===
! Manipulate( type(VectorValue_T) quantity, &
!     type(VectorValue_T) a, type(VectorValue_T) b, rv c, char* mstr, &
!       &  log spreadflag, char*dimList )
! === (end of api) ===
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

  type arrayTemp_T
     real(rv), dimension(:,:), pointer :: VALUES => NULL() ! shaped like a
  end type arrayTemp_T
  
  type(arrayTemp_T), dimension(:), save, pointer :: primitives => null()

  logical, parameter :: COUNTEMPTY = .true.
  logical, parameter :: DEEBUG     = .false.                 ! Usually FALSE
  integer, parameter, public :: MAXMANIPULATIONLEN = 128     ! Is this enough?
  integer, parameter, public :: NO_ERROR_CODE = 0

contains ! =====     Public Procedures     =============================

  subroutine Manipulate( QUANTITY, A, B, C, STR, &
    & SPREADFLAG, DIMLIST )
    ! Args:
    type (VectorValue_T), intent(inout) :: QUANTITY
    type (VectorValue_T), pointer :: A
    type (VectorValue_T), pointer :: B
    real(rv)                      :: C          ! constant "c" in manipulation
    character (len=*), intent(in) :: STR        ! manipulation encoded as a string
    logical, intent(in)           :: SPREADFLAG ! ignore shape, mask, etc.
    character(len=*), intent(in)  :: DIMLIST ! E.g., 's' to shift surfaces, not chans
    ! Evaluate mstr assuming it's of the form
    ! expr1 [op1 expr2]
    ! where each expr is either a primitive 'x' (one of {a, b, or c})
    ! or else ['('] 'x op y' [')']
    ! where 'op' is one of {+, -, *, /,<,>}

    ! Method:
    ! Progressively collapse all the '(..)' pairs into their values
    ! (stored in primitives database)
    ! until only primitives remain
    ! Then evaluate the primitives
    !
    ! Limitations:
    ! Does not check for unmatched parens or other illegal syntax

    ! Improvements to be made:
    ! (1) Check for illegal syntax 
    ! (2) Make ops into array, and loop over them where convenient

    ! Internal variables
    logical :: ALREADYVALPARENS
    character (len=MAXMANIPULATIONLEN) :: mstr ! manipulation being manipulated
    integer, parameter :: MAXNESTINGS=64 ! Max number of '(..)' pairs
    character(len=MAXMANIPULATIONLEN) :: collapsedstr
    integer :: IFUN
    integer :: level
    logical :: MAPFUNCTION
    integer :: np ! number of primitives
    character(len=MAXMANIPULATIONLEN) :: part1
    character(len=MAXMANIPULATIONLEN) :: part2
    character(len=MAXMANIPULATIONLEN) :: part3
    character(len=4) :: vchar
    integer, parameter :: NFUNNAMES = 26
    character(len=8), dimension(NFUNNAMES) :: FUNCOLONS
    character(len=8), dimension(NFUNNAMES), parameter :: FUNNAMES = &
      & (/ 'stddev  ', 'rms     ', 'median  ', 'mean    ', 'max     ', &
      &    'min     ', 'count   ', 'slip    ', 'shift   ', 'channel ', &
      &    'surface ', 'instance', 'height  ', 'lon     ', 'lat     ', &
      &    'sza     ', 'map     ', 'log10   ', 'log     ', 'ln      ', &
      &    'exp     ', 'ifpos   ', 'ifneg   ', 'sign    ', 'abs     ' ,&
      &    'trunc   '/)
      
    ! logical, parameter :: DEEBUG = .true.
    ! Executable
    do iFun = 1, NFUNNAMES
      FUNCOLONS(iFun) = trim(FunNames(iFun)) // ':'
    enddo
    mstr = str
    if ( DeeBUG ) print *, 'mstr: ', trim(mstr)
    MapFunction = ( index(mstr, 'map' ) > 0 )
    nullify(primitives)
    np = 0


    ! Find any terms composed of digits (i.e., literal numbers) ddd and
    ! mark each as val(ddd)
    call markDigits( lowerCase(mstr), collapsedstr )
    if ( DEEBUG ) call outputNamedValue( 'collapsedstr', collapsedstr )

    mstr = collapsedstr
    ! Replace 'e-' with 'e_' to avoid splitting fortran numeric notation
    call ReplaceSubString( mstr, collapsedstr, 'e-', 'e_', &
      & which='all', no_trim=.true. )
    mstr = collapsedstr

    ! We're unable to ensure operator precedence
    ! so we'll attempt to identify multiplications and divisions
    ! and surround such subexpressions with extra parentheses

    call reorderPrecedence( mstr, collapsedstr )
    if ( DeeBUG ) then
      print *, 'incoming ', mstr
      print *, 'after reordering precedence ', collapsedstr
    endif
    mstr = collapsedstr

    ! 1st--make sure spaces surround each operator
    ! (It takes two steps for each to avoid threat of infinite loop)
    call ReplaceSubString( mstr, collapsedstr, '+', ' & ', &
      & which='all', no_trim=.true. )
    call ReplaceSubString( collapsedstr, mstr, '&', '+', &
      & which='all', no_trim=.true. )

    call ReplaceSubString( mstr, collapsedstr, '*', ' & ', &
      & which='all', no_trim=.true. )
    call ReplaceSubString( collapsedstr, mstr, '&', '*', &
      & which='all', no_trim=.true. )

    call ReplaceSubString( mstr, collapsedstr, '-', ' & ', &
      & which='all', no_trim=.true. )
    call ReplaceSubString( collapsedstr, mstr, '&', '-', &
      & which='all', no_trim=.true. )

    call ReplaceSubString( mstr, collapsedstr, '/', ' & ', &
      & which='all', no_trim=.true. )
    call ReplaceSubString( collapsedstr, mstr, '&', '/', &
      & which='all', no_trim=.true. )

    call ReplaceSubString( mstr, collapsedstr, '<', ' & ', &
      & which='all', no_trim=.true. )
    call ReplaceSubString( collapsedstr, mstr, '&', '<', &
      & which='all', no_trim=.true. )

    call ReplaceSubString( mstr, collapsedstr, '>', ' & ', &
      & which='all', no_trim=.true. )
    call ReplaceSubString( collapsedstr, mstr, '&', '>', &
      & which='all', no_trim=.true. )

    call ReplaceSubString( mstr, collapsedstr, '^', ' & ', &
      & which='all', no_trim=.true. )
    call ReplaceSubString( collapsedstr, mstr, '&', '^', &
      & which='all', no_trim=.true. )

    collapsedstr = lowerCase(mstr)
    if ( DEEBUG ) call outputNamedValue( 'collapsedstr', collapsedstr )

    ! Restore 'e-'
    mstr = collapsedstr
    call ReplaceSubString( mstr, collapsedstr, 'e_', 'e-', &
      & which='all', no_trim=.true. )
    alreadyValParens = .false.
    ! Collapse every sub-formula nested within parentheses
    do level =1, MAXNESTINGS ! To prevent endlessly looping if ill-formed
      if ( index( collapsedstr, '(' ) < 1 ) exit
      collapsedstr = stretch( collapsedstr, options='o[(]' )
      call SplitNest ( collapsedstr, part1, part2, part3 )
      ! For Pete's sake, is this a bug you're coding around?
      if ( ncopies( trim(part2), ':' ) > 1 ) then
        call ParensForAll ( collapsedstr, mstr )
        collapsedstr = mstr
        cycle
      elseif ( index( trim(part2), 'val:' ) > 1 .and. .not. alreadyValParens ) then
        alreadyValParens = .true.
        mstr = ParenthesizeVal( collapsedstr )
        if ( DeeBUG ) then
          print *, 'Adding parens to val ', part2
          print *, 'before ', trim(collapsedstr)
          print *, 'after ', trim(mstr)
        endif
        collapsedstr = mstr
        cycle
      elseif ( any( indexes( trim(part2), FunColons) > 1 ) .and. &
        &  .not. alreadyValParens ) then
        alreadyValParens = .true.
        iFun = FindFirst( indexes( trim(part2), FunColons) > 1 )
        mstr = ParenthesizeFun( collapsedstr, trim(FunNames(iFun)) )
        if ( DeeBUG ) then
          print *, 'Adding parens to Fun ', part2, FunNames(iFun)
          print *, 'before ', trim(collapsedstr)
          print *, 'after ', trim(mstr)
        endif
        collapsedstr = mstr
        cycle
      endif
      alreadyValParens = .false.
      ! Now evaluate the part2
      if ( DeeBUG ) then
        print *, 'collapsedstr ', collapsedstr
        print *, 'part1 ', trim(part1)
        print *, 'part2 ', trim(part2)
        print *, 'part3 ', trim(part3)
      endif
      if ( part2 == ' ' ) then
        ! This should never happen with well-formed formulas
        collapsedstr = part1
        cycle
      else
        np = evaluatePrimitive( trim(adjustl(part2)), &
          & a, b, c, &
          & spreadflag, dimList )
        write(vChar, '(i4)') np
      endif
      if ( DeeBUG ) then
        print *, 'vchar ', vchar
      endif
      ! And substitute its value for the spaces it occupied
      if (  part1 // part3 == ' ' ) then
        collapsedstr = vChar
      elseif (  part1 == ' ' ) then
        if ( index(part3, ':') > 1 .and. index(part3, '(') < 1 ) then
          part3 = Parenthesize(part3)
        endif
        ! collapsedstr = trim(vChar) // ' ' // part3
        collapsedstr = catTwoOperands( trim(vChar), part3 )
      elseif ( part3 == ' ' ) then
        ! collapsedstr = trim(part1) // ' ' // vChar
        collapsedstr = catTwoOperands( trim(part1),  vChar )
      else
        ! collapsedstr = trim(part1) // ' ' // trim(vChar) // &
        !   & ' ' // part3
        collapsedstr = catTwoOperands( &
          & trim( catTwoOperands( trim(part1),  trim(vChar) ) ), &
          & part3 )
        if ( DeeBUG ) then
          print *, 'part1 ', trim(part1)
          print *, 'vchar ', trim(part2)
          print *, 'part3 ', trim(part3)
        endif
      endif
      if ( DeeBUG ) then
        print *, 'collapsedstr ', trim(collapsedstr)
      endif
      ! Before exiting, will we need to reorder any priorities?
      if ( index( collapsedstr, '(' ) < 1 ) then
        call reorderPrecedence( collapsedstr, mstr )
        collapsedstr = mstr
        if ( DeeBUG ) then
          print *, 'collapsedstr (reordered) ', trim(collapsedstr)
        endif
      endif
    enddo
    ! Presumably we have collapsed all the nested '(..)' pairs by now
    np = evaluatePrimitive( trim(collapsedstr), &
      & a, b, c, &
      & spreadflag, dimList )
    if ( DeeBUG ) then
      print *, 'np ', np
      print *, 'size(database) ', size(primitives)
    endif
    if ( .not. associated ( quantity%mask ) ) then
      quantity%values = 0.
    else
      where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
        quantity%values = 0.
      end where
    endif
    if ( np < 1 .or. np > size(primitives) ) then
      print *, 'np ', np
      print *, 'size(database) ', size(primitives)
      call Announce_Error ( 'Illegal index for primitives array' )
      return
    endif
    if ( spreadFlag ) then
      ! Ignores mask, shape, etc. if we set the "spread" field
      do level=1, min( size(quantity%values, 1), size(a%values, 1) )
        quantity%values(level, :) = primitives(np)%values(level, 1)
      enddo
      if ( level <=  size(quantity%values, 1) ) &
        quantity%values(level:, :) = primitives(np)%values(level-1, 1)
    elseif ( MapFunction ) then
      ! Ignores mask, shape, etc. if we used the "map" function
      call ReshapeVectorValue( quantity, sourceValues=primitives(np)%values )
    elseif ( .not. associated ( quantity%mask ) ) then
      quantity%values = primitives(np)%values
    else
      where ( iand ( ichar(quantity%mask(:,:)), m_fill ) == 0 )
        quantity%values = primitives(np)%values
      end where
    end if
    if ( DeeBUG ) call dumpPrimitives(primitives)
    call destroyPrimitives(primitives)
  end subroutine Manipulate

  !============ Private procedures ===============
  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( ExtraMessage )
    character (len=*), intent(in), optional :: ExtraMessage

    ! fillerror = max(fillerror,1)
    if ( present(extraMessage) ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & trim(extraMessage) )
    else
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'Calling ANNOUNCE_ERROR' )
    endif
    call output ( " command caused an unrecognized programming error", advance='yes' )
    if ( present(ExtraMessage) )  call output(ExtraMessage, advance='yes')
  end subroutine ANNOUNCE_ERROR

  subroutine doStatFun( qvalue, name, avalues )
    real(rv), intent(out) :: qvalue
    character(len=*), intent(in) :: name ! of the statistical function
    real(rv), dimension(:), intent(in) :: avalues
      select case ( name )
      case ( 'count(a)' )
        qvalue = mlscount( avalues )
      case ( 'min(a)' )
        qvalue = mlsmin( avalues )
      case ( 'max(a)' )
        qvalue = mlsmax( avalues )
      case ( 'mean(a)' )
        qvalue = mlsmean( avalues )
      case ( 'median(a)' )
        qvalue = mlsmedian( avalues )
      case ( 'rms(a)' )
        qvalue = mlsrms( avalues )
      case ( 'stddev(a)' )
        qvalue = mlsstddev( avalues )
      case default
        ! Should not have come here
      end select
  end subroutine doStatFun

  subroutine markDigits( instr, outstr )
    ! Find each instance of a number, composed of consecultive digits
    ! and mark it
    ! E.g., if
    ! instr =  '0.5 * height(a)'
    ! outstr = ' val($0.5) * height(a)'
    ! Args
    character(len=*), intent(in)  :: instr
    character(len=*), intent(out) :: outstr
    ! Internal variables
    ! logical, parameter            :: DEEBUG = .true.
    character(len=1)              :: c
    integer                       :: i          ! char num of instr
    integer                       :: e          ! char num of outstr
    logical                       :: gotDigit
    logical                       :: gotExponent
    character(len=*), parameter :: dlist='1234567890.' ! These are digits
    character(len=*), parameter :: elist='eE'          ! Fortran adds these
    character(len=*), parameter :: signs='+-'          ! legal after 'e'
    ! Executable
    if ( DEEBug ) print *, 'instr ', instr
    outstr = instr
    e = 0
    gotDigit = .false.
    gotExponent = .false.
    do i = 1, len_trim(instr)
      c = instr(i:i)
      if ( index(dlist, c ) > 0 ) then
        ! This was a digit: was it the first?
        if ( gotDigit ) then
          ! Nope, we are just lengthening our number
          e = e + 1
          outstr(e:e) = c
        else
          ! This is the first digit of a number
          ! Distinguish it from index into primitives db
          ! by use of 'val' function and '$' marker
          outstr(e+1:e+5) = 'val($'
          e = e + 6
          outstr(e:e) = c
        endif
        gotDigit = .true.
        gotExponent = .false.
      elseif ( gotExponent ) then
        ! Check that we're not using fortran's '4.9e-6' notation
        if ( index(signs, c ) > 0 ) then
          ! With Fortran notation, we are just lengthening our exponent
          e = e + 1
          outstr(e:e) = c
        else
          ! We have come to the end of our digits
          outstr(e+1:e+1) = ')'
          e = e + 2
          outstr(e:e) = c
          gotDigit = .false.
        endif
      elseif ( gotDigit ) then
        ! Check that we're not using fortran's '4.9e-6' notation
        if ( index(elist, c ) > 0 ) then
          ! With Fortran notation, we are just lengthening our number
          gotExponent = .true.
          e = e + 1
          outstr(e:e) = c
        else
          ! We have come to the end of our digits
          outstr(e+1:e+1) = ')'
          e = e + 2
          outstr(e:e) = c
          gotDigit = .false.
        endif
      else
        e = e + 1
        outstr(e:e) = c
      endif
    enddo
    ! Did we ever close the parentheses?
    if ( gotDigit ) outstr(e+1:e+1) = ')'
    if ( DeeBug ) print *, 'outstr ', outstr
  end subroutine markDigits

  subroutine reorderPrecedence( mstr, collapsedstr )
    ! Identify all the terms where each term are separated by
    ! the lower-precedence operators {+, -,<,>}
    ! If any terms contain higher-precedence operators {*, /}
    ! then surround them by parentheses
    character(len=*), intent(in)  :: mstr
    character(len=*), intent(out) :: collapsedstr
    ! Internal variables
    integer :: i
    integer :: n
    character(len=(len(mstr)+3)) :: element
    character(len=(len(mstr)+3)) :: temp
    ! Executable
    ! Nah--just return
    ! collapsedstr = mstr
    ! return
    ! 1st -- replace each '-' with '+-'
    ! (Don't worry--we'll undo this before returning)
    ! (It takes two steps for each to avoid threat of infinite loop)
    call ReplaceSubString( mstr, collapsedstr, '-', '&', &
      & which='all', no_trim=.true. )
    call ReplaceSubString( collapsedstr, temp, '&', '+-', &
      & which='all', no_trim=.true. )

    call ReplaceSubString( temp, collapsedstr, '<', '&', &
      & which='all', no_trim=.true. )
    call ReplaceSubString( collapsedstr, temp, '&', '+<', &
      & which='all', no_trim=.true. )

    call ReplaceSubString( temp, collapsedstr, '>', '&', &
      & which='all', no_trim=.true. )
    call ReplaceSubString( collapsedstr, temp, '&', '+>', &
      & which='all', no_trim=.true. )
    ! Now loop over terms
    n = NumStringElements( temp, COUNTEMPTY, inseparator='+' )
    if ( n < 1 ) then
      call ReplaceSubString( temp, collapsedstr, '+-', ' -', &
        & which='all', no_trim=.false. )
      call ReplaceSubString( collapsedstr, temp, '+<', ' <', &
        & which='all', no_trim=.false. )
      call ReplaceSubString( temp, collapsedstr, '+>', ' >', &
        & which='all', no_trim=.false. )
      call ReplaceSubString( collapsedstr, temp, '+', ' +', &
        & which='all', no_trim=.false. )
      collapsedstr = temp
      return
    endif
    collapsedstr = ' '
    do i=1, n
      call GetStringElement ( temp, element, i, countEmpty, inseparator='+' )
      ! Surround term with parentheses if it's a product or quotient
      ! but not if it's (already) parenthetical or parentheses are not balanced
      if ( ( index(element, '*') > 0 .or. index(element, '/') > 0 .or. &
        & index(element, '^') > 0 ) .and. &
        & .not. isParenthetical(element) .and. isBalanced(element) ) then
        if ( DEEBUG ) print *, 'element to be parenthesized: ', trim(element)
        element = '(' // trim(element) // ')'
      endif
      collapsedstr = catLists( collapsedstr, element, inseparator='+' )
    enddo

    ! Now undo change by reverting all '+-'
    ! (including any that may have been split by parentheses
    call ReplaceSubString( collapsedstr, temp, '+-', ' -', &
      & which='all', no_trim=.false. )
    call ReplaceSubString( temp, collapsedstr, '+(-', '-(', &
      & which='all', no_trim=.false. )

    call ReplaceSubString( collapsedstr, temp, '+<', '<', &
      & which='all', no_trim=.false. )
    call ReplaceSubString( temp, collapsedstr, '+(<', '<(', &
      & which='all', no_trim=.false. )

    call ReplaceSubString( collapsedstr, temp, '+>', '>', &
      & which='all', no_trim=.false. )
    call ReplaceSubString( temp, collapsedstr, '+(>', '>(', &
      & which='all', no_trim=.false. )
    call ReplaceSubString( collapsedstr, temp, '+', ' +', &
      & which='all', no_trim=.false. )
    collapsedstr = temp
  end subroutine reorderPrecedence

  subroutine destroyPrimitives(primitives)
    ! deallocate all the arrays we created
    type(arrayTemp_T), dimension(:), pointer :: primitives
    integer :: i
    if ( .not. associated(primitives) ) return
    if ( size(primitives) < 1 ) return
    do i=1, size(primitives)
      call deallocate_test( primitives(i)%values, &
        & 'values', ModuleName // '/destroyPrimitives' )
    enddo
  end subroutine destroyPrimitives

  subroutine dumpAPrimitive(primitive)
    ! dump all the values in the array
    type(arrayTemp_T), intent(in) :: primitive
    if ( .not. associated(primitive%values) ) then
      call output( 'values not associated ', advance='yes' )
      return
    endif
    if ( size(primitive%values) < 1 ) then
      call output( 'values array is of 0 size ', advance='yes' )
      return
    endif
    call dump( primitive%values, 'values' )
  end subroutine dumpAPrimitive

  subroutine dumpPrimitives(primitives)
    ! dump all the arrays we created
    type(arrayTemp_T), dimension(:), pointer :: primitives
    integer :: i
    if ( .not. associated(primitives) ) then
      call output( 'database not associated ', advance='yes' )
      return
    endif
    if ( size(primitives) < 1 ) then
      call output( 'empty database ', advance='yes' )
      return
    endif
    call output( 'size of primitives database: ', advance='no' )
    call output( size(primitives), advance='yes' )
    do i=1, size(primitives)
      call output ( ' index of primitive: ', advance='no' )
      call output ( i, advance='yes' )
      call dumpAPrimitive( primitives(i) )
    enddo
  end subroutine dumpPrimitives

  integer function AddPrimitiveToDatabase( DATABASE, ITEM )

    ! This function adds a primitive data type to a database of said types,
    ! creating a new database if it doesn't exist.  The result value is
    ! the size -- where it was put.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    ! Dummy arguments
    type (arrayTemp_T), dimension(:), pointer :: DATABASE
    type (arrayTemp_T), intent(in) :: ITEM

    ! Local variables
    type (arrayTemp_T), dimension(:), pointer :: tempDatabase
    !This include causes real trouble if you are compiling in a different
    !directory.
    include "addItemToDatabase.f9h" 

    AddPrimitiveToDatabase = newSize
  end function AddPrimitiveToDatabase

  function evaluatePrimitive( STR, A, B, C, &
    & SPREADFLAG, DIMLIST ) result(VALUE)
    ! Evaluate an expression composed entirely of
    ! (0) constants ('c')
    ! (1) primitives (e.g., '2')
    ! (2) unary operators ('-')
    ! (3) binary operators {'+', '-', '*', '/','<','>'}
    ! (4) recognized functions {'map:', 'exp:', ..}
    ! Dummy args
    character(len=*)                :: STR
    integer                         :: VALUE
    type (VectorValue_T), pointer   :: A
    type (VectorValue_T), pointer   :: B
    real(rv) :: C                     ! constant "c" in manipulation
    logical, intent(in)             :: SPREADFLAG ! ignore shape, mask, etc.
    character(len=*), intent(in)    :: DIMLIST ! E.g., 's' to shift surfaces, not chans
    ! Internal variables
    logical, dimension(3)           :: average
    character(len=3)                :: avg
    integer, parameter              :: cc = 1      ! channel
    integer, parameter              :: ss = cc + 1 ! surface
    integer, parameter              :: ii = ss + 1 ! instance
    ! logical, parameter              :: DEEBUG = .true.
    logical                         :: done
    ! fun is blank unless a prior one left us "hungry" for an arg
    character(len=8)                :: fun ! {'exp', 'log', etc.}
    integer                         :: elem
    logical                         :: hit
    integer                         :: iChannel
    integer                         :: ind
    integer                         :: instance
    integer                         :: isurf
    integer                         :: k  ! counter for 1st qty index
    character(len=3)                :: lastOp ! {'+', '-', '*', '/'}
    integer                         :: n
    logical                         :: negating
    type (arrayTemp_T)              :: newone
    integer                         :: NoChans
    integer                         :: NoInstances
    integer                         :: NoSurfs
    character(len=8)                :: op
    type (arrayTemp_T)              :: part
    integer                         :: partID
    real(rv)                        :: qvalue
    integer, dimension(2)           :: shp
    integer                         :: surf
    character(len=32)               :: variable
    ! Executable
    average = (/ .true., .false., .false. /)
    if ( len_trim(dimList) > 0 ) then
      average(cc) = index(dimList, 'c') > 0
      average(ss) = index(dimList, 's') > 0
      average(ii) = index(dimList, 'i') > 0
    endif
    avg = ' '
    if ( average(cc) ) avg = trim(avg) // 'c'
    if ( average(ss) ) avg = trim(avg) // 's'
    if ( average(ii) ) avg = trim(avg) // 'i'
    shp = shape(a%values)
    call allocate_test( newone%values, shp(1), shp(2), &
      & 'values', ModuleName // '/evaluatePrimitive(newone)' )
    call allocate_test( part%values, shp(1), shp(2), &
      & 'part%values', ModuleName // '/evaluatePrimitive(part)' )

    if ( deeBug ) then
      print *, 'Complete dump of database'
      call dumpPrimitives(primitives)
    endif

    done = .false.
    negating = .false.
    elem = 0
    lastOp = 'nul' ! 'or'
    newone%values = 0.
    n = NumStringElements( trim(str), countEmpty=.false., &
      & inseparator=' ' )
    if ( DeeBUG ) then
      print *, n, ' str: ', trim(str)
    endif
    partID = -1
    fun = ' '
    hit = .false.
    do
      ! go through the elements, re-evaluating every time we "hit" a primitive
      ! Otherwise revising our lastOp or negating status
      elem = elem + 1
      call GetStringElement ( trim(str), variable, elem, &
        & countEmpty=.false., inseparator=' ' )
      if ( DeeBUG ) then
        print *, elem, ' variable: ', trim(variable)
      endif
      select case( trim(variable) )
      case ('a')
        partID = -1
        part%values = a%values
        hit = .true.
      case ('b')
        partID = -2
        part%values = b%values
        hit = .true.
      case ('c')
        partID = -3
        part%values = c
        hit = .true.
      case ('+')
        lastOp = '+'
        hit = .false.
      case ('*')
        lastOp = '*'
        hit = .false.
      case ('/')
        lastOp = '/'
        hit = .false.
      case ('^')
        lastOp = '^'
        hit = .false.
      case ('-') ! could be unary or binary; how do we tell?
        if ( hit ) then ! already have a primitive; looking for an op
          lastOp = '-'
          hit = .false.
        else
          ! case ('not', '~')
          negating = .true.
          hit = .false.
        endif
      case ('<')
        lastOp = '<'
        hit = .false.
      case ('>')
        lastOp = '>'
        hit = .false.
      case (' ')
        call Announce_Error ( 'parse error of:' // trim(str) )
      case default
        ind = index(variable, ':')
        if ( deeBug ) print *, 'ind of ":" ', ind
        if ( ind > 1 ) then
          ! A function name
          fun = variable(1:ind-1)
          hit = .false.
        elseif ( index(variable, '$') > 0 ) then
          ! A literal number
          if ( deeBug ) print *, 'Trying to read number from ' // variable
          variable = adjustl(variable)
          read( variable(2:), * ) qvalue
          part%values = qvalue
          hit = .true.
          if ( deeBug ) then
            print *, 'part"s values after ' // trim(lastOp) // trim(variable)
            call dumpAPrimitive(part)
          endif
        else
          ! An index into the primitives db
          if ( deeBug ) print *, 'Trying to read partID from ' // variable
          read( variable, * ) partID
          if ( partID < 1 ) then
            print *, 'partID: ', partID
            call Announce_Error ( 'partID too small' )
            return
          elseif( partID > size(primitives) ) then
            print *, 'partID: ', partID
            call Announce_Error ( 'partID too big' )
            return
          endif
          part%values = primitives(partID)%values
          hit = .true.
          if ( deeBug ) then
            print *, 'part"s values after ' // trim(lastOp) // trim(variable)
            call dumpAPrimitive(part)
            print *, 'based on'
            call dumpAPrimitive(primitives(partID))
          endif
        endif
      end select
      if ( hit ) then
        if ( negating ) part%values = -part%values
        op = lastOp
        if ( fun /= ' ' ) op = fun
        if ( deebug ) then
          print *, 'Carry out a ', trim(op)
        endif
        select case(op)
        case ('nul')
            newone%values = part%values
        case ('+')
            newone%values = newone%values + part%values
        case ('-')
            ! print *, 'Subtracting'
            newone%values = newone%values - part%values
        case ('*')
            newone%values = newone%values * part%values
        case ('/')
          where ( part%values /= 0._rv )
            newone%values = newone%values / part%values
          end where
        case ('^')
          where ( newone%values > 0._rv )
            newone%values = newone%values ** part%values
          elsewhere
            newone%values = 0.
          end where
        case ('<')
            newone%values = min( newone%values, part%values )
        case ('>')
            newone%values = max( newone%values, part%values )
        ! Now the functions
        case ('val')
            newone%values = part%values
        case ('abs')
            newone%values = abs( part%values )
        case ('sign')
            where ( part%values /= 0._rv )
              newone%values = sign(1._rv, part%values)
            end where
        case ('trunc')
            newone%values = aint(part%values)
        case ('ifpos')
            where ( part%values > 0._rv )
              newone%values = 1._rv
            end where
        case ('ifneg')
            where ( part%values < 0._rv )
              newone%values = 1._rv
            end where
        case ('exp')
            newone%values = exp( part%values )
        case ('log', 'ln')
            where ( part%values > 0._rv )
              newone%values = log(part%values)
            elsewhere
              newone%values = 0.
            end where
        case ('log10')
            where ( part%values > 0._rv )
              newone%values = log10(part%values)
            elsewhere
              newone%values = 0.
            end where
        ! You can use map to reshape quantities with equal total
        ! size, but distributed differently among channels, sirfs, instances
        ! However, in this context, newone and part are both primitives
        ! and are shaped identically
        case ('map')
            newone%values = part%values
            ! call ReshapeVectorValue( newone, sourceValues=part%values )
            ! call output( 'Calling function map', advance='yes' )
        case ( 'channel', 'surface', 'instance', 'height', 'lon', 'lat', 'sza', &
          & 'phi', 'soltime', 'losangle', 'hid' )
          ! These might be useful for filling arrays with indexes
          NoChans     = a%template%NoChans
          NoInstances = a%template%NoInstances
          NoSurfs     = a%template%NoSurfs
          newone%values = 1
          if ( NoChans*NoSurfs*NoInstances < 2 ) cycle
          do instance=1, NoInstances
            do iSurf=1, NoSurfs
              surf = 1
              if ( .not. a%template%stacked ) surf = iSurf
              do iChannel=1, NoChans
                select case(op)
                case ('channel')
                  newone%values(iChannel + (isurf-1)*NoChans, instance) = &
                    & iChannel
                case ('surface')
                  newone%values(iChannel + (isurf-1)*NoChans, instance) = &
                    & iSurf
                case ('height')
                  if ( a%template%coherent ) then
                    newone%values(iChannel + (isurf-1)*NoChans, instance) = &
                      & a%template%surfs(iSurf, 1)
                  else
                    newone%values(iChannel + (isurf-1)*NoChans, instance) = &
                      & a%template%surfs(iSurf, instance)
                  endif
                case ('instance')
                  newone%values(iChannel + (isurf-1)*NoChans, instance) = &
                    & instance
                case ('lat')
                  newone%values(iChannel + (isurf-1)*NoChans, instance) = &
                    & a%template%GeodLat(surf, instance)
                case ('lon')
                  newone%values(iChannel + (isurf-1)*NoChans, instance) = &
                    & a%template%lon(surf, instance)
                case ('sza')
                  newone%values(iChannel + (isurf-1)*NoChans, instance) = &
                    & a%template%solarZenith(surf, instance)
                case ('soltime')
                  newone%values(iChannel + (isurf-1)*NoChans, instance) = &
                    & a%template%solarTime(surf, instance)
                case ('losangle')
                  newone%values(iChannel + (isurf-1)*NoChans, instance) = &
                    & a%template%losAngle(surf, instance)
                case ('phi')
                  newone%values(iChannel + (isurf-1)*NoChans, instance) = &
                    & a%template%phi(surf, instance)
                case ('hid')
                  newone%values(iChannel + (isurf-1)*NoChans, instance) = &
                    & tai93s2hid(a%template%time(surf, instance))
                end select                
              enddo
            enddo
          enddo
        case ('shift', 'slip')
          ! These are useful for recurrence relations, frequency shifts,
          ! and applying any filter that spans multiple heights or
          ! multiple instances
          ! By default, they apply only to the fastest changing index
          ! in this order: channel, surface, instance
          ! shift(a[c,s,i]) = a[c+1,s,i]
          ! slip(a[c,s,i])  = a[c-1,s,i]
          ! If dimList='s'
          ! shift(a[c,s,i]) = a[c,s+1,i]
          ! If dimList='i'
          ! shift(a[c,s,i]) = a[c,s,i+1]
          
          NoChans     = max(a%template%NoChans, 1)
          NoInstances = a%template%NoInstances
          NoSurfs     = a%template%NoSurfs
          if ( DEEBUG ) then
            call outputnamedValue( 'NoChans      ', NoChans      )
            call outputnamedValue( 'NoInstances  ', NoInstances  )
            call outputnamedValue( 'NoSurfs      ', NoSurfs      )
            call outputnamedValue( 'op           ', op           )
            call outputnamedValue( 'dimList      ', dimList      )
            call outputnamedValue( 'shape(newone)', shape(newone%values)      )
            call outputnamedValue( 'shape(part)  ', shape(part%values)      )
          endif
          if ( NoChans*NoSurfs*NoInstances < 2 ) cycle
          select case(op)
          case ('shift')
            if ( NoChans > 1 .or. index(dimList,'c') > 0 ) then
              do iSurf = 1, NoSurfs
                newone%values(1 + (isurf-1)*NoChans:isurf*NoChans-1, :) = &
                  & part%values(2 + (isurf-1)*NoChans:isurf*NoChans, :)
              enddo
              newone%values(NoChans*NoSurfs, :) = part%values(1, :) ! cyclic
            elseif ( NoSurfs > 1 .or. index(dimList,'s') > 0 ) then
              if ( DEEBUG ) call outputnamedValue( 'shape(newone)  ', &
                & shape(newone%values(1 :1+noChans*(NoSurfs-1)-1:noChans, :)) )
              if ( DEEBUG ) call outputnamedValue( 'shape(part)  ', &
                & shape(part%values(1+noChans :1+noChans*noChans*NoSurfs-1:noChans, :)) )
              do iChannel = 1, NoChans
                newone%values(iChannel :iChannel+noChans*(NoSurfs-1)-1:noChans, :) = &
                  & part%values(iChannel+noChans :iChannel+noChans*NoSurfs-1:noChans, :)
              enddo
              newone%values(NoChans*NoSurfs, :) = part%values(1, :) ! cyclic
            else
              newone%values(:, 1:NoInstances-1) = &
                & part%values(: , 2:NoInstances)
              newone%values(:, NoInstances) = part%values(:, 1) ! cyclic
            endif
            ! call output( 'Calling function shift', advance='yes' )
          case ('slip')
            if ( NoChans > 1 .or. index(dimList,'c') > 0 ) then
              do iSurf = 1, NoSurfs
                newone%values(2 + (isurf-1)*NoChans:isurf*NoChans, :) = &
                  & part%values(1 + (isurf-1)*NoChans:isurf*NoChans-1, :)
              enddo
              newone%values(1, :) = part%values(NoChans*NoSurfs, :) ! cyclic
            elseif ( NoSurfs > 1 .or. index(dimList,'s') > 0 ) then
              do iChannel = 1, NoChans
                newone%values(iChannel+noChans :iChannel+noChans*NoSurfs-1:noChans, :) = &
                  & part%values(iChannel :iChannel+noChans*(NoSurfs-1)-1:noChans, :)
              enddo
              newone%values(1, :) = part%values(NoChans*NoSurfs, :) ! cyclic
            else
              newone%values(:, 2:NoInstances) = &
                & part%values(: , 1:NoInstances-1)
              newone%values(:, 1) = part%values(:, NoInstances) ! cyclic
            endif
            ! call output( 'Calling function slip', advance='yes' )
          end select                
        ! statistical function cases
        case ( 'count', 'min', 'max', 'mean', 'median', 'rms', 'stddev' )
          ! These are harder--we must interpret how to gather
          ! or "average" the data
          ! By default we average over channels, not surfaces or instances
          ! but dimList may cause us to pick out
          ! different choises; e.g., dimList='c,s' tells us to average
          ! over both channels and surfaces
          NoChans     = a%template%NoChans
          NoInstances = a%template%NoInstances
          NoSurfs     = a%template%NoSurfs
          ! if ( dontSumHeights .and. dontSumInstances ) then
          if ( DEEBUG ) then
            call outputNamedValue( 'dimList', dimList )
            call outputNamedValue( 'average', average )
            call outputNamedValue( 'avg', avg )
          endif
          select case (avg)
          case ( 'c' )
          ! if ( .not. average(ss) .and. .not. average(ii) ) then
            ! Average over channels only
            if ( DEEBUG ) call output( 'Average over channels only', advance='yes' )
            do instance = 1, NoInstances
              do iSurf = 1, NoSurfs
                call doStatFun( qvalue, &
                  & trim(op) // '(a)', &
                  & part%values(1+(iSurf-1)*NoChans:iSurf*NoChans, instance) )
                if ( spreadFlag ) then
                  newone%values(1+(iSurf-1)*NoChans:iSurf*NoChans, instance) = qvalue
                else
                  newone%values(1+(iSurf-1)*NoChans, instance) = qvalue
                endif
              enddo
            enddo
          case ( 'cs' )
          ! elseif ( dontSumInstances ) then
          ! elseif ( .not. average(ii) ) then
            if ( DEEBUG ) call output( 'Average over channels and surfaces', advance='yes' )
            do instance = 1, NoInstances
              call doStatFun( qvalue, trim(op) // '(a)', &
                & part%values(:, instance) )
              if ( spreadFlag ) then
                newone%values(:, instance) = qvalue
              else
                newone%values(1, instance) = qvalue
              endif
            enddo
          case ( 'ci' )
          ! elseif ( dontSumHeights ) then
          ! elseif ( .not. average(ss) ) then
            if ( DEEBUG ) call output( 'Average over channels and instances', advance='yes' )
            do iSurf = 1, NoSurfs
              call doStatFun( qvalue, trim(op) // '(a)', &
                & part%values(iSurf, :) )
              if ( spreadFlag ) then
                newone%values(iSurf, :) = qvalue
              else
                newone%values(iSurf, 1) = qvalue
              endif
            enddo
          case ( 'csi' )
          ! else
            ! Average over everything: c, s, and i
            if ( DEEBUG ) call output( 'Average over everything: c, s, and i', advance='yes' )
            select case ( op )
            case ( 'min' )
              qvalue = mlsmin( part%values )
            case ( 'max' )
              qvalue = mlsmax( part%values )
            case ( 'mean' )
              qvalue = mlsmean( part%values )
            case ( 'median' )
              qvalue = mlsmedian( part%values )
            case ( 'rms' )
              qvalue = mlsrms( part%values )
            case ( 'stddev' )
              qvalue = mlsstddev( part%values )
            case default
              ! Should not have come here
            end select
            if ( spreadFlag ) then
              newone%values = qvalue
            else
              newone%values(1, 1) = qvalue
            endif
          case ( 'si' )
            if ( DEEBUG ) call output( 'Average over instances and surfaces', advance='yes' )
            do iChannel = 1, NoChans
              call doStatFun( qvalue, trim(op) // '(a)', &
                & reshape( &
                &   part%values(iChannel:iChannel + noChans*(noSurfs-1):noChans, :), &
                &  (/ NoSurfs*NoInstances /) ) )
              if ( spreadFlag ) then
                newone%values(iChannel:iChannel + noChans*(noSurfs-1):noChans, :) = qvalue
              else
                newone%values(iChannel, :) = qvalue
              endif
            enddo
          ! endif
          case ( 's' )
            if ( DEEBUG ) call output( 'Average over only surfaces', advance='yes' )
            do instance = 1, NoInstances
              do iChannel = 1, NoChans
                call doStatFun( qvalue, trim(op) // '(a)', &
                  & part%values(iChannel:iChannel + noChans*(noSurfs-1):noChans, instance) )
                if ( spreadFlag ) then
                  newone%values(iChannel:iChannel + noChans*(noSurfs-1):noChans, instance) = qvalue
                else
                  newone%values(iChannel, instance) = qvalue
                endif
              enddo
            enddo
          case ( 'i' )
            if ( DEEBUG ) call output( 'Average over only instances', advance='yes' )
            k = 0
            do iSurf = 1, NoSurfs
              do iChannel = 1, NoChans
                k = k + 1
                call doStatFun( qvalue, trim(op) // '(a)', &
                  & part%values(k, :) )
                if ( spreadFlag ) then
                  newone%values(k, 1) = qvalue
                else
                  newone%values(k, :) = qvalue
                endif
              enddo
            enddo
          case default
            ! How could this happen?
          end select
        case default
          ! How could this happen?
            call MLSMessage( MLSMSG_Error, ModuleName, &
              & avg // ' not recognized as a list of averageable dims' )
        end select
        fun = ' '
        negating = .false.
        if ( deeBug ) then
          print *, 'newone"s values after ' // trim(lastOp) // trim(variable)
          call dumpAPrimitive(newone)
        endif
      endif
      if ( DeeBUG ) then
        print *, 'variable ', variable
        print *, 'partID ', partID
        print *, 'hit ', hit
        print *, 'negating ', negating
        print *, 'lastOp ', lastOp
      endif
      done = ( elem >= n )
      if ( done ) exit
    enddo
    value = AddPrimitiveToDatabase( primitives, newone )
    call deallocate_test( part%values, &
      & 'part%values', ModuleName // '/evaluatePrimitive' )
    if ( .not. DEEBUG ) return
    print *, 'value ', value
    print *, 'newone"s values ' // trim(str)
    call dumpAPrimitive(newone)

    print *, 'values stored in db '
    call dumpAPrimitive(primitives(value))
  end function evaluatePrimitive

  function catTwoOperands( part1, part2 ) result ( str )
    ! cat together part1 and part2 with a space between them
    ! unless the last non-blank character of part1
    ! and the 1st non-blank character of part2 aren't operators
    ! in which case put a colon ':' between them
    ! E.g., if 
    ! part1 = 'a + b'
    ! and part2 = '/ c' then str = 'a + b / c'
    ! but if part1 = 'map'
    ! and part2 = 'c - a' then str = 'map: c - a'
    ! args
    character(len=*), intent(in)           :: part1
    character(len=*), intent(in)           :: part2
    character(len=MAXMANIPULATIONLEN)        :: str
    ! internal variables
    character(len=1), dimension(9), parameter :: ops = &
      &          (/ '+', '-', '*', '/' , '(', ')', '^', '<', '>' /)
    character(len=1), parameter :: oparen='(', cparen=')'
    ! character(len=1), parameter :: oparen=' ', cparen=' '
    character(len=1) :: part1Tail, part2Head
    integer :: maxind
    integer :: n
    integer, dimension(4) :: inds
    ! logical, parameter :: DEEBUG = .true.
    ! Executable
    n = max(1, len_trim(part1))
    part1Tail = part1(n:n)
    part2Head = adjustl(part2)
    if ( all( indexes(part1Tail // part2Head, ops) == 0 ) ) then
      ! Mark function name "hungry" for an arg by adding a ':' to
      ! str = trim(part1) // ': ' // adjustl(part2)
      ! Must also check if part 1 contains an embedded operator
      ! Or an embedded '('
      inds = indexes( part1, (/ '+', '-', '*', '/' /), mode='last' )
      maxind = maxval(inds)
      if ( maxind < 1 ) then
        if ( DeeBUG ) then
          print *, 'cat with part1 ', trim(part1)
          print *, 'cat with part2 ', trim(part2)
        endif
        if ( index(part1, '(') > 0 ) then
          str = trim(part1) // ': ' // adjustl(part2)
          return
        endif
        str = oparen // trim(part1) // ': ' // trim(adjustl(part2) ) // cparen
      elseif ( index(part1, '(') > 0 ) then
        str = trim(part1) // ': ' // adjustl(part2)
        return
      else
        str = part1(:maxind) // oparen // trim(part1(maxind+1:)) &
          & // ': ' // trim(adjustl(part2) ) // cparen
      endif
    else
      str = trim(part1) // ' ' // adjustl(part2)
    endif
  end function catTwoOperands

  function isParenthetical ( str ) result ( itIs )
    ! TRUE if 1st non-blank is '(' and last non-blank is ')'
    ! Do we really need to make this a recursive function?
    character(len=*), intent(in) :: str
    logical                      :: itIs
    itIs = index( adjustl(str), '('           ) == 1 .and. &
      &    index( trim(str), ')', back=.true. ) == len_trim(str)
  end function isParenthetical

  function Parenthesize ( str ) result ( newStr )
    ! Add parentheses to a function
    character(len=*), intent(in) :: str
    character(len=len(str)+8)    :: newStr
    ! Method: starting from the ':', go backward to the first blank
    ! and forward to the first blank
    ! put '( ' or ') ' in place of them
    ! Internal variables
    integer :: col
    integer :: dbInd
    integer :: firstBl
    integer :: lastBl
    newStr = str
    col = index( str, ':' )
    if ( col < 1 ) return
    firstBl = findLast( str(:col), ' ' )
    dbInd = col+1
    if ( str(col+1:col+1) == ' ' ) dbInd = col+2
    lastBl = FindFirst( str(dbInd:), ' ' )
    if ( firstBl > 1 ) then
      newStr = str(:firstBl-1) // '( ' // str(firstBl+1:dbInd+lastBl-1) // ' )'
    else
      newStr = '( ' // str(firstBl+1:dbInd+lastBl-1) // ' )'
    endif
    if ( dbInd+lastBl-1 < len_trim(str) ) &
      newStr = trim(newStr) // str(dbInd+lastBl:)
  end function Parenthesize

  function ParenthesizeFun ( str, Fun ) result ( newStr )
    ! Add parentheses to a Fun function
    ! E.g., turn ' ... Fun: 4 ...'
    ! to ' ... (Fun: 4) ...'
    character(len=*), intent(in) :: str
    character(len=*), intent(in) :: Fun
    character(len=len(str)+8)    :: newStr
    ! Method: starting from the 'fun:', go forward to the first blank
    ! put '( ' or ') ' in place of them
    ! Internal variables
    integer :: col
    integer :: dbInd
    integer :: firstBl
    integer :: lastBl
    integer :: lenFunName
    ! Executable
    lenFunName = len(Fun)
    newStr = str
    col = index( str, Fun // ':' )
    if ( col < 1 ) return
    firstBl = findFirst( str(col+lenFunName+2:), ' ' )
    if ( col > 1 ) then
      newStr = str(:col-1) // '( ' // Fun // ':'
    else
      newStr =                '( ' // Fun // ':'
    endif
    newStr = trim(newStr) // ' ' // str(col+lenFunName+2:col+lenFunName+2+firstBl-1) // &
      & ')' // str(col+lenFunName+2+firstBl:)
  end function ParenthesizeFun

  function ParenthesizeVal ( str ) result ( newStr )
    ! Add parentheses to a val function
    ! E.g., turn ' ... val: 4 ...'
    ! to ' ... (val: 4) ...'
    character(len=*), intent(in) :: str
    character(len=len(str)+8)    :: newStr
    ! Method: starting from the 'val:', go forward to the first blank
    ! put '( ' or ') ' in place of them
    ! Internal variables
    integer :: col
    integer :: dbInd
    integer :: firstBl
    integer :: lastBl
    newStr = str
    col = index( str, 'val:' )
    if ( col < 1 ) return
    firstBl = findFirst( str(col+5:), ' ' )
    if ( col > 1 ) then
      newStr = str(:col-1) // '( ' // 'val:'
    else
      newStr =                '( ' // 'val:'
    endif
    newStr = trim(newStr) // ' ' // str(col+5:col+5+firstBl-1) // &
      & ')' // str(col+5+firstBl:)
  end function ParenthesizeVal

  subroutine ParensForAll ( str, outStr )
    ! Args
    character(len=*), intent(in)  :: str
    character(len=*), intent(out) :: outstr
    ! Internal variables
    integer :: i
    integer :: ia, ib
    integer :: n    ! number of ':'-separated elements
    integer :: NLOCBI
    character(len=1) :: a, b
    character(len=len(str)+8), dimension(24) :: stringArray
    character(len=len(str)+8)                :: stringElm
    ! logical, parameter                       :: DEEBUG = .true.
    ! Executable
    n = NumStringElements( str, countEmpty, &
      & ':' )
    if ( n < 2 ) return
    call list2Array( str, stringArray, countEmpty, ':', &
      & IgnoreLeadingSpaces = .false. )
    do i=1, n-1
      a = ' '
      b = ' '
      ! Find loc of last blank of element i
      nlocBi = FindLast( trim(stringArray(i)), ' ' )
      if ( nLocBi < 2 ) cycle
      ! a is last non-blank char before nLocBi? Is a '('?
      ia = FindLast( trim(stringArray(i)(:nLocBi)), ' ', reverse=.true. )
      if ( ia > 0 ) a = stringArray(i)(ia:ia)
      ! b is second non-blank char after start of element i+1? Is a ')'?
      ib = FindNext( trim(stringArray(i+1)), ' ', current=1, reverse=.true. )
      if ( ib > 0 ) then
        b = stringArray(i+1)(ib:ib)
      else
        ib = len_trim(stringArray(i+1))
      endif
      if ( DEEBUG ) then
        print *, 'stringArray(i) ', trim(stringArray(i))
        print *, 'stringArray(i+1) ', trim(stringArray(i+1))
        print *, 'ia ', ia
        print *, 'ib ', ib
      endif
      if ( a == '(' .or. b == ')' ) cycle
      if ( ia > 0 ) then
        stringArray(i) = stringArray(i)(:ia) // ' (' // stringArray(i)(ia+1:)
      else
        stringArray(i) = ' (' // stringArray(i)
      endif
      stringArray(i+1) = stringArray(i+1)(:ib) // ' )' // stringArray(i+1)(ib+1:)
    enddo
    call Array2List( stringArray(1:n), outStr,  inseparator=':' )
  end subroutine ParensForAll

  function isBalanced ( str ) result ( itIs )
    use MLSStrings, only: ncopies
    ! Are numbers of '(' and ')' in str equal?
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
    ! If 0 when done, we;'re balanced, not if not
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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------
end module ManipulationUtils
!=============================================================================

!
! $Log$
! Revision 2.17  2015/05/05 17:47:04  pwagner
! MAXMANIPULATIONLEN made public
!
! Revision 2.16  2015/04/21 17:49:38  pwagner
! Added trunc() function
!
! Revision 2.15  2015/03/28 02:48:36  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.14  2014/09/05 01:11:15  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.13  2014/09/05 00:49:07  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 2.12  2014/07/18 23:17:32  pwagner
! Aimed for consistency in names passed to allocate_test
!
! Revision 2.11  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.10  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.9  2013/08/06 23:02:24  pwagner
! Needed another change to handle nested functions
!
! Revision 2.8  2013/07/13 00:17:08  pwagner
! Fixed another bug; this one affected numbers
!
! Revision 2.7  2013/06/17 21:59:13  pwagner
! Fixed serious bugs
!
! Revision 2.6  2013/04/22 17:51:01  pwagner
! Added soltime, phi, losangle, hid to Gelocation functions
!
! Revision 2.5  2013/01/18 01:43:48  pwagner
! Fixed various bugs impacting nested expressions
!
! Revision 2.4  2013/01/14 21:21:06  pwagner
! Fixed certain bugs; others may require major changes
!
! Revision 2.3  2012/11/14 01:00:04  pwagner
! Use dimList for choosing which of {csi} to average over
!
! Revision 2.2  2012/11/08 23:33:48  pwagner
! dimList field lets us specifiy whether to shift by [c,s,i]
!
! Revision 2.1  2012/10/09 00:49:12  pwagner
! First commit
!
