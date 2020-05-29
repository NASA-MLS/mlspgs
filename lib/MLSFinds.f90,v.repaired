! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MLSFinds
  use MLSStrings_0, only: LowerCase

! These are just the simplest Finds

  implicit none
  private

! === (start of toc) ===                                                 
!     c o n t e n t s                                                    
!     - - - - - - - -                                                    

!     (subroutines and functions)
! FindAll       Find all logicals in the array that are true, or all the
!               integers in the array equal to the probe, "matches"
! FindFirst     Find the first logical in the array that is true, or the
!               first [integer,real,double] in the array equal to the probe
! FindLast      Find the last instead
! FindLongestRange
!               Find the longest stretch of consecutive matches
! FindUnique    Return only the unique elements of a set;
!                 formally, a set contains only unique elements, so this
!                 will in fact reduce any improper set to a proper set
! FindNext      Find the next instead
! FindLongestStretch    Find the longest stretch of integers
! === (end of toc) ===

! === (start of api) ===
! FindAllCharacter (char* set(:), char* it, int which(:), [int how_many], &
!                  [char* re_mainder(:)], [int which_not(:)])      
! FindAllInteger (int set(:), int it, int which(:), [int how_many], &
!                  [int re_mainder(:)], [int which_not(:)])      
! FindAllLogical (log set(:), int which(:), [int how_many], &
!                  [int which_not(:)])      
! int FindFirstCharacter (char* set(:), char* probe[, char* Tol] )      
! int FindFirstInCharacter ( char* Str, char* set(:), int c1 )      
! int FindFirstNumType (numtype set(:), numtype probe, [numtype tol], &
!      [numtype Period])
!     (where numtype can be an int, real, or dble)
! int FindFirstLogical (log condition(:))      
! int FindFirstSubString (char* set, char* probe, [log reverse])      
! int FindFirstVector (int set(:,:), int probe(:), [log reverse])      
! int FindFirstRun (int set(:), int probe(:), [log reverse])      
! int FindLastCharacter (char* set(:), char* probe, [log reverse])
! int FindLastNumType (numtype set(:), numtype probe, [numtype tol], &
!      [numtype Period], [log reverse])
!     (where numtype can be an int, real, or dble)
! int FindLastLogical (log condition(:), [log reverse])    
! int FindLastSubString (char* set, char* probe, [log reverse])
! FindLongestCharacter (char* set(:), char* probe, int range(2))
! FindLongestInteger (int set(:), int probe, int range(2))
! FindLongestLogical (log condition(:), int range(2))
! FindLongestSubString (char* set, char* probe, int range(2), [log reverse])
! FindLongestStretch( int which(:), int how_many, int range(2) )
! FindUniqueCharacter (char* set(:), char* unique(:), [int nUnique], &
!    [int(:) counts])
! FindUniqueCharacterSubString (char* set, char* unique, [int nUnique], &
!    [int(:) counts])
! FindUniqueInteger (int set(:), int unique(:), [int nUnique], [int(:) counts])
! The optional argument reverse, which may be supplied to some of the functions
! causes the usual matching condition to be reversed; e.g. FindLast
! with reverse=TRUE returns the last index of the array that does not
! match the probe

! === (end of api) ===

  public :: Findall, FindFirst, FindLast, FindLongestRange, &
    & Findnext, FindUnique, FindLongeststretch, N_Indices

  interface FindFirst
    module procedure FindFirstInteger, FindFirstlogical
    module procedure FindFirstReal, FindFirstdouble, findFirstlogical2d
    module procedure FindFirstSubstring, findFirstcharacter, findFirstIncharacter
    module procedure FindFirstVector
    module procedure FindFirstRun
  end interface

  interface FindLast
    module procedure FindLastInteger, FindLastlogical, FindLastcharacter
    module procedure FindLastReal, FindLastdouble, FindLastlogical2d
    module procedure FindLastSubstring
  end interface

  interface FindLongestRange
    module procedure FindLongestInteger, FindLongestlogical, FindLongestcharacter
    module procedure FindLongestSubstring
  end interface

  interface FindUnique
    module procedure FindUniqueInteger, FindUniquecharacter
    module procedure FindUniqueCharactersubstring
    module procedure FindUniqueReal, FindUniquedouble
  end interface

  interface FindNext
    module procedure FindNextInteger, FindNextlogical, FindNextcharacter
    module procedure FindNextSubstring
  end interface

  interface FindAll
    module procedure FindAllInteger, FindAlllogical, FindAllcharacter
    module procedure FindAllSubstring
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  integer, parameter                   :: MaxListLength = 500

contains ! =====     Public Procedures     =============================

  ! ---------------------------------------------  FindAll  -----
  ! This family of subroutines Finds all the matching elements
  ! not just the first, next, or last
  ! Return which i of set[i] = it
  ! If character-valued,
  !   (case-sensitive, ignores trailing blanks, but alert to leading blanks)
  ! Optionally, may return also:
  ! how_many    of them do
  ! which_not   of them (which don't)
  ! re_mainder  of them != it
  ! e.g. given set = /(4, 3, 1, 2, 1, 3 )/ and it = 1
  ! produces which = /(3, 5)/, 
  !      which_not = /(1, 2, 4, 6)/, 
  !       how_many = 2,
  !     re_mainder = /(4, 3, 2, 3)/

  ! Note that which and which_not are arrays of indices
  ! while re_mainder is an array of elements

  ! This may be useful,
  ! e.g. in dump for reshaping an array to suppress any dims 
  ! that are identically 1

  ! Maybe this should be rewritten to use fraternal functions FindFirst
  ! or FindNext

  ! ---------------------------------------------  FindAllCharacter  -----
  subroutine FindAllCharacter ( set, it, which, how_many, re_mainder, which_not )
    ! Formal arguments
    character(len=*), intent(in), dimension(:)   :: set
    character(len=*), intent(in)                 :: it
    integer, intent(out), dimension(:)           :: which
    integer, intent(out), optional               :: how_many
    character(len=*), intent(out), dimension(:), optional :: re_mainder
    integer, intent(out), dimension(:), optional :: which_not

    ! local variables
    integer :: i, i_which, i_re_mainder
    
    if ( size(set) < 1 .or. size(which) < 1 ) then
      if ( present(how_many) ) how_many = 0
      if ( present(re_mainder) ) re_mainder = ''
      return
    end if
    i_which = 0
    i_re_mainder = 0
    do i=1, size(set)
      if ( trim(set(i)) == trim(it) ) then
        i_which = i_which+1
        which(min(size(which), i_which)) = i
      else
        i_re_mainder = i_re_mainder+1
        if ( present(which_not) ) &
          & which_not(min(size(which_not), i_re_mainder)) = i
        if ( present(re_mainder) ) &
          & re_mainder(min(size(re_mainder), i_re_mainder)) = set(i)
      end if
    end do
    if ( present(how_many) ) how_many = i_which

  end subroutine FindAllCharacter

  ! ---------------------------------------------  FindAllInteger  -----
  subroutine FindAllInteger ( set, it, which, how_many, re_mainder, which_not )
    ! Formal arguments
    integer, intent(in), dimension(:)  ::           set
    integer, intent(in)                ::           it
    integer, intent(out), dimension(:) ::           which
    integer, intent(out), optional ::               how_many
    integer, intent(out), dimension(:), optional :: re_mainder
    integer, intent(out), dimension(:), optional :: which_not

    ! local variables
    integer :: i, i_which, i_re_mainder
    
    if ( size(set) < 1 .or. size(which) < 1 ) then
      if ( present(how_many) ) how_many = 0
      if ( present(re_mainder) ) re_mainder = 0
      return
    end if
    i_which = 0
    i_re_mainder = 0
    do i=1, size(set)
      if ( set(i) == it ) then
        i_which = i_which+1
        which(min(size(which), i_which)) = i
      else
        i_re_mainder = i_re_mainder+1
        if ( present(which_not) ) &
          & which_not(min(size(which_not), i_re_mainder)) = i
        if ( present(re_mainder) ) &
          & re_mainder(min(size(re_mainder), i_re_mainder)) = set(i)
      end if
    end do
    if ( present(how_many) ) how_many = i_which

  end subroutine FindAllInteger

  ! ---------------------------------------------  FindAllLogical  -----
  subroutine FindAllLogical ( set, which, how_many, which_not )
    ! Return which i of set[i] = .true.
    ! optionally, may return also:
    ! how many of them do
    ! which_not of them (which don't)
    ! e.g. given set = /(F, F, T, F, T, F )/
    ! produces which = /(3, 5)/, 
    !      which_not = /(1, 2, 4, 6)/, 
    !       how_many = 2,
    
    ! Note that which and which_not are arrays of array indices
    ! Formal arguments
    logical, intent(in), dimension(:)  ::           set
    integer, intent(out), dimension(:) ::           which
    integer, intent(out), optional ::               how_many
    integer, intent(out), dimension(:), optional :: which_not

    ! local variables
    integer :: i, i_which, i_re_mainder
    
    if ( size(set) < 1 .or. size(which) < 1 ) then
      if ( present(how_many) ) how_many = 0
      return
    end if
    i_which = 0
    i_re_mainder = 0
    do i=1, size(set)
      if ( set(i) ) then
        i_which = i_which+1
        which(min(size(which), i_which)) = i
      else
        i_re_mainder = i_re_mainder+1
        if ( present(which_not) ) &
          & which_not(min(size(which_not), i_re_mainder)) = i
      end if
    end do
    if ( present(how_many) ) how_many = i_which

  end subroutine FindAllLogical

  ! ---------------------------------------------  FindAllSubString  -----
  subroutine FindAllSubString ( set, it, which, how_many, re_mainder, which_not )
    ! Return which substring i of set[i:i] = it
    ! Formal arguments
    character(len=*), intent(in)                 :: set
    character(len=*), intent(in)                 :: it
    integer, intent(out), dimension(:)           :: which
    integer, intent(out), optional               :: how_many
    character(len=*), intent(out), optional      :: re_mainder
    integer, intent(out), dimension(:), optional :: which_not
    character(len=*), parameter                  :: method = 'new'
    if ( method == 'new' ) then
      call newFindAllSubString ( set, it, which, how_many, re_mainder, which_not )
    else
      call oldFindAllSubString ( set, it, which, how_many, re_mainder, which_not )
    endif

  end subroutine FindAllSubString

  ! -------------------------------------------  FindFirstCharacter  -----
  integer function FindFirstCharacter ( Set, Probe, Tol )
    ! Find the first element in the array Set that matches Probe
    ! By default
    ! case-sensitive, ignores trailing blanks, but alert to leading blanks
    ! These can be overridden by Tol string
    !           If Tol
    !   contains      effect
    !   --------      ------
    !      c          case-insensitive
    !      p          partial match OK, if probe contained entirely in  set
    !      s          partial match OK, if set entirely contained in  probe
    character(len=*), dimension(:), intent(in) :: Set
    character(len=*), intent(in) :: Probe
    character(len=*), optional, intent(in) :: Tol ! 
    ! Internal variables
    character(len=8)          :: options
    character(len=len(probe)) :: p
    character(len=len(set))   :: test

    ! Executable code
    options = ' '
    if ( present(Tol) ) options = Tol
    do FindFirstCharacter = 1, size(set)
      if ( index( options, 'c' ) > 0 ) then
        test = lowercase(set(FindFirstCharacter))
        p = lowercase(probe)
      else
        test = set(FindFirstCharacter)
        p = probe
      endif
      if ( index( options, 'p' ) > 0 ) then
        if ( index( test, trim(p)) > 0 ) return
      elseif ( index( options, 's' ) > 0 ) then
        if ( index( p, trim(test)) > 0 ) return
      elseif ( trim(test) == trim(p) ) then
        return
      endif
    end do
    FindFirstCharacter = 0
  end function FindFirstCharacter

  ! -------------------------------------------  FindFirstInCharacter  -----
  integer function FindFirstInCharacter ( Str, Set, c1 )
    ! Find the first occurrence in str of any element in the array Set
    character(len=*), intent(in)               :: Str
    character(len=*), dimension(:), intent(in) :: Set
    integer, intent(out)                       :: c1 ! its substring index
    ! Internal variables
    integer                   :: j
    integer                   :: k

    ! Executable code
    c1 = len_trim(str) + 1
    FindFirstInCharacter = 0
    do k = 1, size(set)
      j = index( str, trim(set(k)) )
      if ( j > 0 ) then
        ! Does this elem occur before the earliest we found?
        if ( j < c1 ) then
          c1 = j
          FindFirstInCharacter = k
        endif
      endif
    end do
  end function FindFirstInCharacter

  ! ----------------------  FindFirst[Integer,real,double]  -----
  function FindFirstInteger ( Set, Probe, Tol, Period ) Result( theFirst )
    ! Find the first element in the array Set that is equal to Probe
    integer, dimension(:), intent(in) :: Set
    integer, intent(in)               :: Probe
    integer, intent(in), optional     :: Tol  ! Ignored; purely for consistency
    integer, intent(in), optional     :: Period  ! 
    integer                           :: theFirst

    ! Executable code
    if ( present(Period) ) then
      ! If we're checking for coincidences allowing for a periodic function
      ! then we want to allow the differences to be close within an
      ! integer multiple of the period
      do theFirst = 1, size(set)
        if ( mod( set(theFirst) - probe, period ) == 0 ) return
      end do
      theFirst = 0
      return
    end if

    do theFirst = 1, size(set)
      if ( set(theFirst) == probe ) return
    end do

    theFirst = 0
  end function FindFirstInteger

  function FindFirstReal ( Set, Probe, Tol, Period ) Result( theFirst )
    ! Find the first element in the array Set that is equal to Probe
    ! Or whose difference < Tol
    integer, parameter :: RK = kind(1.0e0)
    real(rk), dimension(:), intent(in) :: Set
    real(rk), intent(in)           :: Probe
    real(rk), intent(in), optional :: Tol
    real(rk), intent(in), optional :: Period

    integer                        :: theFirst
    real(rk)                       :: myTol
    real(rk)                       :: q

    ! Executable code
    myTol = 0.0
    if ( present(tol) ) myTol = abs(tol)
    if ( present(Period) ) then
      ! If we're checking for coincidences allowing for a periodic function
      ! then we want to allow the differences to be close within an
      ! integer multiple of the period
      do theFirst = 1, size(set)
        q = abs( set(theFirst) - probe ) / Period
        if ( q - int(q) <= myTol/Period ) return
      end do
      theFirst = 0
      return
    end if

    do theFirst = 1, size(set)
      if ( abs( set(theFirst) - probe) <= myTol ) return
    end do

    theFirst = 0
  end function FindFirstReal

  function FindFirstDouble ( Set, Probe, Tol, Period ) Result( theFirst )
    ! Find the first element in the array Set that is equal to Probe
    ! Or whose difference < Tol
    integer, parameter :: RK = kind(1.0d0)
    real(rk), dimension(:), intent(in) :: Set
    real(rk), intent(in)           :: Probe
    real(rk), intent(in), optional :: Tol
    real(rk), intent(in), optional :: Period

    integer                        :: theFirst
    real(rk)                       :: myTol
    real(rk)                       :: q

    ! Executable code
    myTol = 0.0
    if ( present(tol) ) myTol = abs(tol)
    if ( present(Period) ) then
      ! If we're checking for coincidences allowing for a periodic function
      ! then we want to allow the differences to be close within an
      ! integer multiple of the period
      do theFirst = 1, size(set)
        q = abs( set(theFirst) - probe ) / Period
        if ( q - int(q) <= myTol/Period ) return
      end do
      theFirst = 0
      return
    end if

    do theFirst = 1, size(set)
      if ( abs( set(theFirst) - probe) <= myTol ) return
    end do

    theFirst = 0
  end function FindFirstDouble

  ! -------------------------------------------  FindFirstLogical  -----
  integer function FindFirstLogical ( condition )
    ! Find the first logical in the array that is true
    logical, dimension(:), intent(in) :: CONDITION

    ! Executable code
    do FindFirstLogical = 1, size(condition)
      if ( condition(FindFirstLogical) ) return
    end do
    FindFirstLogical = 0
  end function FindFirstLogical

  ! -------------------------------------------  FindFirstLogical2D  -----
  integer function FindFirstLogical2D ( condition, options )
    ! Find the first logical in index 1 of the array that is true
    ! for all (any) of the 2nd indexes
    logical, dimension(:,:), intent(in)           :: CONDITION
    character(len=*), optional, intent(in)        :: options
    ! options is a character-valued parameter in the form
    ! '-[t][n][l]'    with the following effect
    !
    !      t   transpose 1st and 2nd indices of array
    !      n   choose sense 'any' while testing condition
    !      l   choose sense 'all' while testing condition (default)
    character(len=3) :: myOptions

    ! Executable code
    myOptions = '-l'
    if ( present(options) ) then
      if ( index(options, 'n') > 0 ) myOptions = '-n'
      if ( index(options, 't') > 0 ) myOptions = myOptions(:2) // 't'
    end if
    select case(myOptions)
    case ( '-l' )
      do FindFirstLogical2D = 1, size(condition, 1)
        if ( all(condition(FindFirstLogical2D, :)) ) return
      end do
    case ( '-n' )
      do FindFirstLogical2D = 1, size(condition, 1)
        if ( any(condition(FindFirstLogical2D, :)) ) return
      end do
    case ( '-lt' )
      do FindFirstLogical2D = 1, size(condition, 2)
        if ( all(condition(:, FindFirstLogical2D)) ) return
      end do
    case ( '-nt' )
      do FindFirstLogical2D = 1, size(condition, 2)
        if ( any(condition(:, FindFirstLogical2D)) ) return
      end do
    case default
      ! This should never occur--treat it as if no matching condition found
    end select
    FindFirstLogical2D = 0
  end function FindFirstLogical2D

  ! -------------------------------------------  FindFirstSubString  -----
  integer function FindFirstSubString ( Set, Probe, reverse, Tol )
    ! Find the first sub-string in the string Set that is (not) equal to Probe
    ! Along with FindLastSubstring, does no more than intrinsic index function
    ! except for optional arg reverse which allows us to return index of
    ! substring (actually only a single char) that does *not* match
    character(len=*), intent(in) :: Set
    character(len=1), intent(in) :: Probe
    logical, optional, intent(in) :: reverse
    character(len=1), optional, intent(in) :: Tol ! ignored; generic consistency
    ! Internal variables
    logical :: myReverse

    ! Executable code
    myReverse = .false.
    if ( present(reverse) ) myReverse = reverse
    do FindFirstSubString = 1, len(set)
      if ( myReverse ) then
        if ( set(FindFirstSubString:FindFirstSubString) /= probe ) return
      else
        if ( set(FindFirstSubString:FindFirstSubString) == probe ) return
      end if
    end do
    FindFirstSubString = 0
  end function FindFirstSubString

  ! -------------------------------------------  FindFirstVector  -----
  integer function FindFirstVector ( Set, Probe, reverse, Tol )
    ! Find the first vector in the Set that is (not) equal to Probe
    ! Note that probe is a 1d array, Set a 2d array, and that
    ! we match each of the first indexes in Set
    ! See also FindFirstRun
    integer, dimension(:,:), intent(in) :: Set
    integer, dimension(:), intent(in)   :: Probe
    logical, optional, intent(in)       :: reverse
    character(len=1), optional, intent(in) :: Tol ! ignored; generic consistency
    ! Internal variables
    logical :: myReverse

    ! Executable code
    myReverse = .false.
    if ( present(reverse) ) myReverse = reverse
    do FindFirstVector = 1, size(Set, 2)
      if ( myReverse ) then
        if ( any(set(:,FindFirstVector) /= probe) ) return
      else
        if ( all(set(:,FindFirstVector) == probe) ) return
      end if
    end do
    FindFirstVector = 0
  end function FindFirstVector

  ! -------------------------------------------  FindFirstRun  -----
  integer function FindFirstRun ( Set, Probe, reverse, Tol )
    ! Find the first Run in the Set that is (not) equal to Probe
    ! Note that probe is a 1d array, Set also a 1d array, and that
    ! we match a run of consecutive indexes in Set
    ! See also FindFirstVector
    integer, dimension(:), intent(in)   :: Set
    integer, dimension(:), intent(in)   :: Probe
    logical, optional, intent(in)       :: reverse
    character(len=1), optional, intent(in) :: Tol ! ignored; generic consistency
    ! Internal variables
    logical :: myReverse
    integer :: n
    ! Executable code
    myReverse = .false.
    if ( present(reverse) ) myReverse = reverse
    n = size(Probe)
    FindFirstRun = 0
    if ( n > size(Set) ) return
    FindFirstRun = 1
    if ( n == size(Set) ) then
      if ( myReverse ) then
        if ( any(set /= probe) ) return
      else
        if ( all(set == probe) ) return
      end if
      FindFirstRun = 0
      return
    endif
    do FindFirstRun = 1, size(Set) - n + 1
      if ( myReverse ) then
        if ( any(set(FindFirstRun:FindFirstRun+n-1) /= probe) ) return
      else
        if ( all(set(FindFirstRun:FindFirstRun+n-1) == probe) ) return
      end if
    end do
    FindFirstRun = 0
  end function FindFirstRun

  ! These next could be done by reversing the list order and
  ! calling findFirst
  ! -------------------------------------------  FindLastCharacter  -----
  function FindLastCharacter ( Set, Probe, Reverse ) result(LAST)
    ! Find the last element in the array Set that is equal to Probe
    ! (case-sensitive, ignores trailing blanks, but alert to leading blanks)
    character(len=*), dimension(:), intent(in) :: Set
    character(len=*), intent(in)               :: Probe
    logical, optional, intent(in)              :: REVERSE
    integer                                    :: LAST
    logical                                    :: myReverse
    ! Executable code
    myReverse = .false.
    if ( present(reverse) ) myReverse = reverse
    if ( myReverse ) then
      do Last = size(set), 1, -1
        if ( trim(set(Last)) /= trim(probe) ) return
      end do
    else
      do Last = size(set), 1, -1
        if ( trim(set(Last)) == trim(probe) ) return
      end do
    endif
    Last = 0
  end function FindLastCharacter

  ! -------------------------------------------  FindLastInteger  -----
  function FindLastInteger ( Set, Probe, Reverse ) result(LAST)
    ! Find the last element in the array Set that is equal to Probe
    integer, dimension(:), intent(in) :: Set
    integer, intent(in) :: Probe
    logical, optional, intent(in)              :: REVERSE
    integer                                    :: LAST
    logical                                    :: myReverse
    ! Executable code
    myReverse = .false.
    if ( present(reverse) ) myReverse = reverse
    if ( myReverse ) then
      do Last = size(set), 1, -1
        if ( set(Last) /= probe ) return
      end do
    else
      do Last = size(set), 1, -1
        if ( set(Last) == probe ) return
      end do
    endif
    Last = 0
  end function FindLastInteger

  function FindLastReal ( Set, Probe, Tol, Reverse ) Result( theLast )
    ! Find the Last element in the array Set that is equal to Probe
    ! Or whose difference < Tol
    integer, parameter :: RK = kind(1.0e0)
    real(rk), dimension(:), intent(in) :: Set
    real(rk), intent(in)           :: Probe
    real(rk), intent(in), optional :: Tol
    logical, intent(in),  optional :: Reverse

    integer                        :: theLast
    real(rk)                       :: myTol
    ! Internal variables
    logical :: myReverse

    ! Executable code
    myReverse = .false.
    if ( present(reverse) ) myReverse = reverse
    myTol = 0.0
    if ( present(tol) ) myTol = abs(tol)

    if ( myReverse ) then
      do theLast = size(set), 1, -1
        if ( abs( set(theLast) - probe) > myTol ) return
      end do
    else
      do theLast = size(set), 1, -1
        if ( abs( set(theLast) - probe) <= myTol ) return
      end do
    endif
    theLast = 0
  end function FindLastReal

  function FindLastDouble ( Set, Probe, Tol, Reverse ) Result( theLast )
    ! Find the Last element in the array Set that is equal to Probe
    ! Or whose difference < Tol
    integer, parameter :: RK = kind(1.0d0)
    real(rk), dimension(:), intent(in) :: Set
    real(rk), intent(in)           :: Probe
    real(rk), intent(in), optional :: Tol
    logical, intent(in),  optional :: Reverse

    integer                        :: theLast
    real(rk)                       :: myTol
    ! Internal variables
    logical :: myReverse

    ! Executable code
    myReverse = .false.
    if ( present(reverse) ) myReverse = reverse
    myTol = 0.0
    if ( present(tol) ) myTol = abs(tol)
    if ( myReverse ) then
      do theLast = size(set), 1, -1
        if ( abs( set(theLast) - probe) > myTol ) return
      end do
    else
      do theLast = size(set), 1, -1
        if ( abs( set(theLast) - probe) <= myTol ) return
      end do
    endif
    theLast = 0
  end function FindLastDouble

  ! -------------------------------------------  FindLastLogical  -----
  function FindLastLogical ( condition, Reverse ) result(LAST)
    ! Find the last logical in the array that is true
    logical, dimension(:), intent(in) :: CONDITION
    logical, optional, intent(in)              :: REVERSE
    integer                                    :: LAST
    logical                                    :: myReverse
    ! Executable code
    myReverse = .false.
    if ( present(reverse) ) myReverse = reverse
    if ( myReverse ) then
      do Last = size(condition), 1, -1
        if ( .not. condition(Last) ) return
      end do
    else
      do Last = size(condition), 1, -1
        if ( condition(Last) ) return
      end do
    endif
    Last = 0
  end function FindLastLogical

  ! -------------------------------------------  FindLastLogical2D  -----
  integer function FindLastLogical2D ( condition, options )
    ! Find the first logical in index 1 of the array that is true
    ! for all (any) of the 2nd indexes
    logical, dimension(:,:), intent(in)           :: CONDITION
    character(len=*), optional, intent(in)        :: options
    ! options is a character-valued parameter in the form
    ! '-[t][n][l]'    with the following effect
    !
    !      t   transpose 1st and 2nd indices of array
    !      n   choose sense 'any' while testing condition
    !      l   choose sense 'all' while testing condition (default)
    character(len=3) :: myOptions

    ! Executable code
    myOptions = '-l'
    if ( present(options) ) then
      if ( index(options, 'n') > 0 ) myOptions = '-n'
      if ( index(options, 't') > 0 ) myOptions = myOptions(:2) // 't'
    end if
    select case(myOptions)
    case ( '-l' )
      do FindLastLogical2D = size(condition, 1), 1 -1
        if ( all(condition(FindLastLogical2D, :)) ) return
      end do
    case ( '-n' )
      do FindLastLogical2D = size(condition, 1), 1, -1
        if ( any(condition(FindLastLogical2D, :)) ) return
      end do
    case ( '-lt' )
      do FindLastLogical2D = size(condition, 2), 1, -1
        if ( all(condition(:, FindLastLogical2D)) ) return
      end do
    case ( '-nt' )
      do FindLastLogical2D = size(condition, 2), 1, -1
        if ( any(condition(:, FindLastLogical2D)) ) return
      end do
    case default
      ! This should never occur--treat it as if no matching condition found
    end select
    FindLastLogical2D = 0
  end function FindLastLogical2D

  ! -------------------------------------------  FindLastSubString  -----
  integer function FindLastSubString ( Set, Probe, reverse )
    ! Find the last substing in the string Set that is (not) equal to Probe
    character(len=*), intent(in) :: Set
    character(len=*), intent(in) :: Probe
    logical, optional, intent(in) :: reverse
    ! Internal variables
    logical :: myReverse

    ! Executable code
    myReverse = .false.
    if ( present(reverse) ) myReverse = reverse
    do FindLastSubString = len(set), 1, -1
      if ( myReverse ) then
        if ( set(FindLastSubString:FindLastSubString) /= probe ) return
      else
        if ( set(FindLastSubString:FindLastSubString) == probe ) return
      end if
    end do
    FindLastSubString = 0
  end function FindLastSubString

  ! -------------------------------------------  FindLongest  -----
  ! Find the Longest stretch of consecutive elements in the array Set
  ! (equal to Probe)
  subroutine FindLongestCharacter ( Set, Probe, Range )
    ! (case-sensitive, ignores trailing blanks, but alert to leading blanks)
    character(len=*), dimension(:), intent(in) :: Set
    character(len=*), intent(in) :: Probe
    integer, dimension(2), intent(out) :: range
    ! Internal variables
    integer, dimension(size(set)) :: which
    integer :: how_many

    ! Executable code
    range = 0
    call FindAll( set, probe, which, how_many )
    call FindLongestStretch( which, how_many, range )
  end subroutine FindLongestCharacter

  subroutine FindLongestInteger ( Set, Probe, Range )
    ! Find the Longest stretch of consecutive elements in the array Set
    ! equal to Probe
    ! (case-sensitive, ignores trailing blanks, but alert to leading blanks)
    integer, dimension(:), intent(in) :: Set
    integer, intent(in) :: Probe
    integer, dimension(2), intent(out) :: range
    ! Internal variables
    integer, dimension(size(set)) :: which
    integer :: how_many

    ! Executable code
    range = 0
    call FindAll( set, probe, which, how_many )
    call FindLongestStretch( which, how_many, range )
  end subroutine FindLongestInteger

  subroutine FindLongestLogical ( Set, Probe, Range )
    ! Find the Longest stretch of consecutive elements in the array Set
    ! equal to Probe
    ! (case-sensitive, ignores trailing blanks, but alert to leading blanks)
    logical, dimension(:), intent(in) :: Set
    logical, intent(in) :: Probe
    integer, dimension(2), intent(out) :: range
    ! Internal variables
    integer, dimension(size(set)) :: which
    integer, dimension(size(set)) :: which_not
    integer :: how_many

    ! Executable code
    range = 0
    call FindAll( set, which, how_many, which_not )
    if ( probe ) then
      call FindLongestStretch( which, how_many, range )
    else
      call FindLongestStretch( which_not, how_many, range )
    end if
  end subroutine FindLongestLogical

  subroutine FindLongestSubString ( Set, Probe, Range, reverse )
    ! Find the Longest stretch of consecutive elements in the array Set
    ! equal to Probe
    ! (case-sensitive, ignores trailing blanks, but alert to leading blanks)
    character(len=*), intent(in)                 :: set
    character(len=1), intent(in)                 :: probe
    integer, dimension(2), intent(out) :: range
    logical, optional, intent(in) :: reverse
    ! Internal variables
    integer, dimension(len(set)) :: which
    integer, dimension(len(set)) :: which_not
    integer :: how_many
    logical :: myReverse

    ! Executable code
    myReverse = .false.
    if ( present(reverse) ) myReverse = reverse
    range = 0
    call FindAll( set, probe, which, how_many, which_not=which_not )
    if ( myReverse ) then
      how_many = len(set) - how_many
      which = which_not
    end if
    call FindLongestStretch( which, how_many, range )
  end subroutine FindLongestSubString
  
   subroutine FindLongestStretch( which, how_many, range )
    ! Given an array of integers, return the longest
    ! stretch of consecutive ones
    integer, dimension(:), intent(in)  :: which
    integer, intent(in)                :: how_many
    integer, dimension(2), intent(out) :: range
    ! Internal variables
    logical, parameter                 :: DEEBug = .false.
    integer                            :: firstWhich
    integer                            :: i
    logical                            :: inSequence
    integer                            :: lastWhich
    integer                            :: stretch
    ! Executable
    range = 0
    if ( how_many < 1 ) return
    range = which(1)
    if ( how_many < 2 ) return
    stretch = 1
    firstWhich = which(1)
    lastWhich = which(1)
    do i=2, how_many
      ! Are we still in sequence?
      if ( DEEBug ) then
        print *, "i is"
        print *, i
        print *, "which(i)"
        print *, which(i)
        print *, "lastWhich"
        print *, lastWhich
      endif
      if ( (which(i)-lastWhich) > 1 ) then
        ! Nope--our sequence ended with lastWhich
        ! How does it compare with the previous record-holder?
        inSequence = .false.
        if ( lastWhich - firstWhich + 1 > stretch ) then
          ! It is longer so it becomes our best bet
          range(1) = firstWhich
          range(2) = lastWhich
          stretch = lastWhich - firstWhich + 1
        ! else
          ! Not longer, so we forget about it
        end if
        firstWhich = which(i)
        lastWhich = which(i)
      else
        ! Yes, we have extended our sequence
        lastWhich = which(i)
        inSequence = .true.
      end if
    end do
    ! Did we end the loop still in sequence?
    if ( inSequence ) then
      if ( lastWhich - firstWhich + 1 > stretch ) then
        ! It is longer so it becomes our new winner
        range(1) = firstWhich
        range(2) = lastWhich
      end if
    end if
  end subroutine FindLongestStretch

  ! --------------------------------------------  FindNextCharacter  -----
  integer function FindNextCharacter ( Set, Probe, Current, &
    & Wrap, Repeat, Reverse )
    ! Find the next element in the array Set that is equal to Probe after the
    ! current one
    ! (case-sensitive, ignores trailing blanks, but alert to leading blanks)
    ! May optionally wrap or repeat
    ! e.g., if wrap is true and current is last true, return first true
    ! e.g., if repeat is true and current is also last true, return current
    ! wrap takes priority over repeat if both are present and true
    character(len=*), dimension(:), intent(in) :: Set
    character(len=*), intent(in) :: Probe
    integer, intent(in) :: Current
    logical, optional, intent(in) :: Wrap
    logical, optional, intent(in) :: Repeat
    logical, optional, intent(in) :: reverse

    ! Local variables
    integer :: I                        ! Loop counter
    logical :: myWrap
    logical :: myRepeat
    logical :: myReverse

    ! Executable code
    myWrap = .false.
    if ( present(wrap) ) myWrap = wrap
    myRepeat = .false.
    if ( present(repeat) ) myRepeat = repeat
    myReverse = .false.
    if ( present(reverse) ) myReverse = reverse
    FindNextCharacter = 0
    ! We'll assume you gave us valid args; otherwise return 0
    if ( current < 1 .or. current > size(set)) return
    ! if ( set(current) /= probe ) return
    ! Now check for current already at end of array
    if ( current < size(set) ) then
      do i = current+1, size(set)
        if ( myReverse ) then
          if ( trim(set(i)) /= trim(probe) ) then
            FindNextCharacter = i
            return
          end if
        else
          if ( trim(set(i)) == trim(probe) ) then
            FindNextCharacter = i
            return
          end if
        end if
      end do
    end if
    ! Uh-oh, this means current is last true
    if ( myWrap ) then
      FindNextCharacter = FindFirst(set,probe)
    else if ( myRepeat ) then
      FindNextCharacter = current
    end if
  end function FindNextCharacter

  ! --------------------------------------------  FindNextInteger  -----
  integer function FindNextInteger ( Set, Probe, Current, Wrap, Repeat )
    ! Find the next element in the array Set that is equal to Probe after the
    ! current one
    ! May optionally wrap or repeat
    ! e.g., if wrap is true and current is last true, return first true
    ! e.g., if repeat is true and current is also last true, return current
    ! wrap takes priority over repeat if both are present and true
    integer, dimension(:), intent(in) :: Set
    integer, intent(in) :: Probe
    integer, intent(in) :: Current
    logical, optional, intent(in) :: Wrap
    logical, optional, intent(in) :: Repeat

    ! Local variables
    integer :: I                        ! Loop counter
    logical :: myWrap
    logical :: myRepeat

    ! Executable code
    myWrap = .false.
    if ( present(wrap) ) myWrap = wrap
    myRepeat = .false.
    if ( present(repeat) ) myRepeat = repeat
    FindNextInteger = 0
    ! We'll assume you gave us valid args; otherwise return 0
    if ( current < 1 .or. current > size(set)) return
    ! if ( set(current) /= probe ) return
    ! Now check for current already at end of array
    if ( current < size(set) ) then
      do i = current+1, size(set)
        if ( set(i) == probe ) then
          FindNextInteger = i
          return
        end if
      end do
    end if
    ! Uh-oh, this means current is last true
    if ( myWrap ) then
      FindNextInteger = FindFirst(set,probe)
    else if ( myRepeat ) then
      FindNextInteger = current
    end if
  end function FindNextInteger

  ! --------------------------------------------  FindNextLogical  -----
  integer function FindNextLogical ( condition, current, wrap, repeat )
    ! Find the next logical in the array that is true after the current one
    ! May optionally wrap or repeat
    ! e.g., if wrap is true and current is last true, return first true
    ! e.g., if repeat is true and current is also last true, return current
    ! wrap takes priority over repeat if both are present and true
    logical, dimension(:), intent(in) :: CONDITION
    integer, intent(in) :: CURRENT
    logical, optional, intent(in) :: WRAP
    logical, optional, intent(in) :: REPEAT

    ! Local variables
    integer :: I                        ! Loop counter
    logical :: myWrap
    logical :: myRepeat

    ! Executable code
    myWrap = .false.
    if ( present(wrap) ) myWrap = wrap
    myRepeat = .false.
    if ( present(repeat) ) myRepeat = repeat
    FindNextLogical = 0
    ! We'll assume you gave us valid args; otherwise return 0
    if ( current < 1 .or. current > size(condition)) return
    ! if ( .not. condition(current)) return
    ! Now check for current already at end of array
    if ( current < size(condition) ) then
      do i = current+1, size(condition)
        if ( condition(i) ) then
          FindNextLogical = i
          return
        end if
      end do
    end if
    ! Uh-oh, this means current is last true
    if ( myWrap ) then
      FindNextLogical = FindFirst(condition)
    else if ( myRepeat ) then
      FindNextLogical = current
    end if
  end function FindNextLogical

  ! --------------------------------------------  FindNextSubString  -----
  integer function FindNextSubString ( Set, Probe, Current, &
    & Wrap, Repeat, Reverse )
    ! Find the next substring in the string Set that is equal to Probe after the
    ! current one
    ! May optionally wrap or repeat
    ! e.g., if wrap is true and current is last true, return first true
    ! e.g., if repeat is true and current is also last true, return current
    ! wrap takes priority over repeat if both are present and true
    character(len=*), intent(in) :: Set
    character(len=1), intent(in) :: Probe
    integer, intent(in) :: Current
    logical, optional, intent(in) :: Wrap
    logical, optional, intent(in) :: Repeat
    logical, optional, intent(in) :: reverse

    ! Local variables
    integer :: I                        ! Loop counter
    logical :: myWrap
    logical :: myRepeat
    logical :: myReverse

    ! Executable code
    myWrap = .false.
    if ( present(wrap) ) myWrap = wrap
    myRepeat = .false.
    if ( present(repeat) ) myRepeat = repeat
    myReverse = .false.
    if ( present(reverse) ) myReverse = reverse
    FindNextSubString = 0
    ! We'll assume you gave us valid args; otherwise return 0
    if ( current < 1 .or. current > len(set)) return
    ! if ( (set(current:current) == probe) .eqv. myReverse ) return
    ! Now check for current already at end of array
    if ( current < len(set) ) then
      do i = current+1, len(set)
        if ( myReverse ) then
          if ( set(i:i) /= probe ) then
            FindNextSubString = i
            return
          end if
        else
          if ( set(i:i) == probe ) then
            FindNextSubString = i
            return
          end if
        end if
      end do
    end if
    ! Uh-oh, this means current is last true
    if ( myWrap ) then
      FindNextSubString = FindFirst(set,probe,reverse)
    else if ( myRepeat ) then
      FindNextSubString = current
    end if
  end function FindNextSubString

  ! --------------------------------------------  FindUnique ------------------
  subroutine FindUniqueInteger ( Set, Unique, nUnique, counts )
    ! Find all unique elements in the array Set
    integer, dimension(:), intent(in)            :: Set
    integer, dimension(:), intent(out)           :: Unique  ! array of unique elements
    integer, optional, intent(out)               :: nUnique ! num of unique elements
    integer, dimension(:), optional, intent(out) :: counts  ! how often each appears

    ! Local variables
    integer :: myUnique(size(Set))
    include 'findunique.f9h'
    Unique = 0
    Unique(1:num) = myUnique(1:num)
    if ( present(counts) ) counts(1:num) = myCounts(1:num)
  end subroutine FindUniqueInteger

  subroutine FindUniqueReal ( Set, Unique, nUnique, counts )
    ! Find all unique elements in the array Set
    real, dimension(:), intent(in)            :: Set
    real, dimension(:), intent(out)           :: Unique  ! array of unique elements
    integer, optional, intent(out)               :: nUnique ! num of unique elements
    integer, dimension(:), optional, intent(out) :: counts  ! how often each appears

    ! Local variables
    real :: myUnique(size(Set))
    include 'findunique.f9h'
    Unique = 0.
    Unique(1:num) = myUnique(1:num)
    if ( present(counts) ) counts(1:num) = myCounts(1:num)
  end subroutine FindUniqueReal

  subroutine FindUniqueDouble ( Set, Unique, nUnique, counts )
    ! Find all unique elements in the array Set
    double precision, dimension(:), intent(in)            :: Set
    double precision, dimension(:), intent(out)           :: Unique  ! array of unique elements
    integer, optional, intent(out)               :: nUnique ! num of unique elements
    integer, dimension(:), optional, intent(out) :: counts  ! how often each appears

    ! Local variables
    double precision :: myUnique(size(Set))
    include 'findunique.f9h'
    Unique = 0.d0
    Unique(1:num) = myUnique(1:num)
    if ( present(counts) ) counts(1:num) = myCounts(1:num)
  end subroutine FindUniqueDouble

  subroutine FindUniqueCharacter ( Set, Unique, nUnique, counts )
    ! Find all unique elements in the array Set
    character(len=*), dimension(:), intent(in)  :: Set
    character(len=*), dimension(:), intent(out) :: Unique   ! array of unique elements
    integer, optional,          intent(out)     :: nUnique  ! num of unique elements
    integer, dimension(:), optional, intent(out) :: counts  ! how often each appears

    ! Local variables
    character(len=len(Set)) :: myUnique(size(Set))
    include 'findunique.f9h'
    Unique = ' '
    Unique(1:num) = myUnique(1:num)
    if ( present(counts) ) counts(1:num) = myCounts(1:num)
  end subroutine FindUniqueCharacter

  subroutine FindUniqueCharacterSubString ( Set, Unique, nUnique, counts )
    ! Find all unique substrings in the character Set
    character(len=*), intent(in)  :: Set
    character(len=*), intent(out) :: Unique   ! composed of unique subvstrings
    integer, optional,          intent(out)     :: nUnique  ! num of unique elements
    integer, dimension(:), optional, intent(out) :: counts  ! how often each appears

    ! Local variables
    integer :: i                        ! Loop counter
    integer :: myCounts(len(Set))
    character(len=len(Set)) :: myUnique
    integer :: num
    integer :: prev

    ! Executable code
    num = 0
    myCounts = 0
    if ( present(nUnique) ) nUnique = 0
    if ( len(Unique) < 1 ) return
    if ( len(Set) < 1 ) then
      return
    elseif ( len(Set) < 2 ) then
      Unique = Set
      if ( present(nUnique) ) nUnique = 1
      if ( present(counts) ) counts = 1
      return
    end if
    num = 1
    myUnique(1:1) = Set(1:1)
    myCounts(1) = 1
    do i=2, len(Set)
      prev = FindFirst( myUnique(1:num), Set(i:i) )
      if ( prev > 0 ) then
        myCounts(prev) = myCounts(prev) + 1
      else
        num = num + 1
        myUnique(num:num) = Set(i:i)
        myCounts(num) = 1
      end if
    end do
    if ( present(nUnique) ) nUnique = num
    num = min( num, len(Unique) )
    Unique = ' '
    Unique(1:num) = myUnique(1:num)
    if ( present(counts) ) counts(1:num) = myCounts(1:num)
  end subroutine FindUniqueCharacterSubString

  ! Private procedures and functions
!============================ Private ==============================
  function n_indices ( str, substr, which ) result ( n )
    ! how many instances of substr occur in str
    character(len=*), intent(in)       :: str
    character(len=*), intent(in)       :: substr
    integer                            :: n
    integer, dimension(:)              :: which
    integer                            :: k
    n = 0
    if ( len_trim(str) == 0 .or. len_trim(substr) == 0 ) return
    ! Method:
    ! We pick off the last instance of substr in whatever remains
    ! of str after picking off the previously last instance.
    ! We do this until we can't find another instance or
    ! whatever remains is blank.
    k = len_trim(str)
    do
      k = index( str(1:k), trim(substr), back=.true. )
      if ( k < 1 ) exit
      n = n + 1
      which(n) = k
      k = k - 1
    enddo
    if ( n < 2 ) return
    ! Reverse the order of indices
    which(1:n) = which(n:1:-1)
  end function n_indices

  subroutine newFindAllSubString ( set, it, which, how_many, re_mainder, which_not )
    ! Return which substring i of set[i:i] = it
    ! Formal arguments
    character(len=*), intent(in)                 :: set
    character(len=*), intent(in)                 :: it
    integer, intent(out), dimension(:)           :: which
    integer, intent(out), optional               :: how_many
    character(len=*), intent(out), optional      :: re_mainder
    integer, intent(out), dimension(:), optional :: which_not
    ! Internal variables
    character(len=len(set))                      :: myremainder
    integer, dimension(MAXLISTLENGTH)            :: myWhich_Not
    integer                                      :: i
    integer                                      :: k
    integer                                      :: pLast
    integer                                      :: m
    integer                                      :: n
    integer                                      :: p
    ! Executable
    m = len_trim(set)
    n = n_indices( set, it, which )
    ! Process any optional args
    if ( present( how_many ) ) how_many = n
    if ( .not. present(re_mainder) .and. .not. present(which_not) ) return
    myRemainder = set
    myWhich_Not(1:m) = (/(k, k=1,m)/)
    if ( n < 1 ) then
      if ( present(which_not) ) which_not = myWhich_not
      if ( present(re_mainder) ) re_mainder  = myremainder
      return
    endif
    k = 0
    plast = 0
    myremainder = ' '
    myWhich_Not = 0
    do i=1, n
      do p=plast+1, which(i)-1
        k = k + 1
        myWhich_not(k) = p
        myremainder(k:k) = set(p:p)
      enddo
      plast = which(i)
    enddo
    if ( which(n) < m ) then
      do p=plast+1, m
        k = k + 1
        myWhich_not(k) = p
        myremainder(k:k) = set(p:p)
      enddo
    endif
    if ( present(which_not) ) which_not = myWhich_not
    if ( present(re_mainder) ) re_mainder  = myremainder
  end subroutine newFindAllSubString 

  subroutine oldFindAllSubString ( set, it, which, how_many, re_mainder, which_not )
    ! Return which substring i of set[i:i] = it
    ! Formal arguments
    character(len=*), intent(in)                 :: set
    character(len=*), intent(in)                 :: it
    integer, intent(out), dimension(:)           :: which
    integer, intent(out), optional               :: how_many
    character(len=*), intent(out), optional :: re_mainder
    integer, intent(out), dimension(:), optional :: which_not

    ! local variables
    integer :: i, i_which, i_re_mainder
    integer :: itl     ! len(it)
    integer :: ip      ! last index starting at i to fit it
    ! Executable
    if ( len(set) < 1 .or. size(which) < 1 ) then
      if ( present(how_many) ) how_many = 0
      if ( present(re_mainder) ) re_mainder = ''
      return
    end if
    itl = len(it)
    i_which = 0
    i_re_mainder = 0
    do i=1, len(set)
      ip = i + itl - 1
      ip = min( ip, len(set) )
      if ( set(i:ip) == it ) then
        i_which = i_which+1
        which(min(size(which), i_which)) = i
      else
        i_re_mainder = i_re_mainder+1
        if ( present(which_not) ) &
          & which_not(min(size(which_not), i_re_mainder)) = i
        if ( present(re_mainder) ) &
          & re_mainder(i_re_mainder:i_re_mainder) = set(i:i)
      end if
    end do
    if ( present(how_many) ) how_many = i_which

  end subroutine oldFindAllSubString

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSFinds

! $Log$
! Revision 2.8  2020/05/29 21:51:13  pwagner
! Implemented speedier version of FindAllSubString
!
! Revision 2.7  2019/11/22 00:36:39  pwagner
! Some timid housekeeping
!
! Revision 2.6  2019/10/16 20:50:53  pwagner
! Added FindFirstInCharacter
!
! Revision 2.5  2019/04/09 20:33:04  pwagner
! Now USEs LowerCase from MLSStrings_0
!
! Revision 2.4  2017/11/03 21:12:06  uid9452
! Must update lastWhich in FindLongestStretch
!
! Revision 2.3  2017/02/11 00:29:42  pwagner
! FindLongestStretch now public
!
! Revision 2.2  2014/07/31 20:20:23  pwagner
! FindFirstCharacter can now find partial match
!
! Revision 2.1  2013/08/12 23:44:05  pwagner
! First commit
!
