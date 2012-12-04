! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MLSSets

! Various operations on sets

! Note:
! A set in mathematics consists of a set of elements
! Each element in a set appears exactly once, i.e. each is unique

! We represent sets with arrays
! We usually do not check whether the elements are in fact unique
! (you can use IsProperSet to check or FindUnique to ensure this) 
! so some operations may not work properly if you supply arrays with
! duplicate elements

! A point about null sets:
! Not all operations have been checked thoroughly for correctness
! when one or more arguments are null sets, i.e. 0-sized
! IsProperSet, IsProperSubset, and IsSubset do give correct results:
! The null set is 
!  (1) a proper set, 
!  (2) a proper subset of any set except itself
!  (3) a subset of itself, but not a proper subset of itself

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
! FindIntersection  
!               Compute indices of elements in intersection of two sets
! FindLast      Find the last instead
! FindLongestRange
!               Find the longest stretch of consecutive matches
! FindNext      Find the next instead
! FindUnique    Return only the unique elements of a set;
!                 formally, a set contains only unique elements, so this
!                 will in fact reduce any improper set to a proper set
! Intersect     Return true if two sets represented by arrays of integers have
!               a common element
! Intersection  Compute intersection of two sets
! IsProperSet   Check that each element is unique
! IsProperSubSet 
!               Check that each element in A is within larger B
! IsSubSet      Check that each element in A is within B
! RelativeComplement
!               Compute complement of set a relative to set b
!               i.e., all elements in b except those in a
! Union         Compute union of two sets
! UnionSize     Compute the size a union of two sets would have.
! === (end of toc) ===

! === (start of api) ===
! FindAllCharacter (char* set(:), char* it, int which(:), [int how_many], &
!                  [char* re_mainder(:)], [int which_not(:)])      
! FindAllInteger (int set(:), int it, int which(:), [int how_many], &
!                  [int re_mainder(:)], [int which_not(:)])      
! FindAllLogical (log set(:), int which(:), [int how_many], &
!                  [int which_not(:)])      
! int FindFirstCharacter (char* set(:), char* probe)      
! int FindFirstNumType (numtype set(:), numtype probe, [numtype tol], &
!      [numtype Period])
!     (where numtype can be an int, real, or dble)
! int FindFirstLogical (log condition(:))      
! int FindFirstSubString (char* set, char* probe, [log reverse])      
! FindIntersection ( set1(:), set2(:), int which1(:), int which2(:),
!      [int how_many] )
! int FindLastCharacter (char* set(:), char* probe, [log reverse])
! int FindLastInteger (int set(:), int probe, [log reverse])
! int FindLastLogical (log condition(:), [log reverse])    
! int FindLastSubString (char* set, char* probe, [log reverse])
! FindLongestCharacter (char* set(:), char* probe, int range(2))
! FindLongestInteger (int set(:), int probe, int range(2))
! FindLongestLogical (log condition(:), int range(2))
! FindLongestSubString (char* set, char* probe, int range(2), [log reverse])
! int FindNextCharacter (char* set(:), char* probe, int current, [log wrap], [log repeat])
! int FindNextInteger (int set(:), int probe, int current, [log wrap], [log repeat]) 
! int FindNextLogical (log condition(:), int current, [log wrap], [log repeat])
! int FindNextSubString (char* set, char* probe, [log wrap], [log repeat], [log reverse])
! FindUniqueCharacter (char* set(:), char* unique(:), [int nUnique], [int(:) counts])
! FindUniqueCharacterSubString (char* set, char* unique, [int nUnique], [int(:) counts])
! FindUniqueInteger (int set(:), int unique(:), [int nUnique], [int(:) counts])
! log Intersect ( set A, set B )
! set Intersection ( set A, set B, [char* options] )
! log isProperSet ( set A )
! log isProperSubset ( set A, set B )
! log isSubset ( set A, set B )
! set RelativeComplement ( set A, set B )
! set Union ( set A, set B )
! int UnionSize ( set A, set B )

! Notes:
! A set is represented either
! as an array of integers, floats, or an array of characters; e.g.
!   logical Intersect ( set a(:), set b(:) )
!   set *Intersection ( set a(:), set b(:) )
!   logical IsProperSet ( set a(:) )
!   logical Is[Proper]SubSet ( set a(:), set b(:) )
!   set *RelativeComplement ( set a(:), set b(:) )
!   set *Union ( set a(:), set b(:) )
!   int UnionSize ( set a(:), set b(:) )

! The optional argument reverse, which may be supplied to some of the functions
! causes the usual matching condition to be reversed; e.g. FindLast
! with reverse=TRUE returns the last index of the array that does not
! match the probe

! === (end of api) ===

  public :: FINDALL, FINDFIRST, FINDINTERSECTION, FINDLAST, FINDLONGESTRANGE, &
    & FINDNEXT, FINDUNIQUE, &
    & INTERSECT, INTERSECTION, ISPROPERSET, ISPROPERSUBSET, ISSUBSET, &
    & RELATIVECOMPLEMENT, &
    & UNION, UNIONSIZE

  interface FindFirst
    module procedure FindFirstInteger, FindFirstLogical, FindFirstCharacter
    module procedure FindFirstReal, FindFirstDouble, FindFirstLogical2D
    module procedure FindFirstSubString
  end interface

  interface FindIntersection
    module procedure FindIntersectionInteger, FindIntersectionReal, &
      & FindIntersectionDouble, FindIntersectionCharacter
  end interface

  interface FindLast
    module procedure FindLastInteger, FindLastLogical, FindLastCharacter
    module procedure FindLastLogical2D
    module procedure FindLastSubString
  end interface

  interface FindLongestRange
    module procedure FindLongestInteger, FindLongestLogical, FindLongestCharacter
    module procedure FindLongestSubString
  end interface

  interface FindNext
    module procedure FindNextInteger, FindNextLogical, FindNextCharacter
    module procedure FindNextSubString
  end interface

  interface FindAll
    module procedure FindAllInteger, FindAllLogical, FindAllCharacter
    module procedure FindAllSubString
  end interface

  interface FindUnique
    module procedure FindUniqueInteger, FindUniqueCharacter
    module procedure FindUniqueCharacterSubString
    module procedure FindUniqueReal, FindUniqueDouble
  end interface

  interface Intersect
    module procedure IntersectInteger, IntersectCharacter
  end interface

  interface Intersection
    module procedure IntersectionInteger, IntersectionCharacter
    module procedure IntersectionReal, IntersectionDouble
  end interface

  interface IsProperSet
    module procedure IsProperSetInteger, IsProperSetCharacter
  end interface

  interface IsProperSubset
    module procedure IsProperSubsetInteger, IsProperSubsetCharacter
  end interface

  interface IsSubSet
    module procedure IsSubSetInteger, IsSubSetCharacter
  end interface

  interface RelativeComplement
    module procedure RelativeComplementInteger, RelativeComplementCharacter
    module procedure RelativeComplementReal, RelativeComplementDouble
  end interface

  interface Union
    module procedure UnionInteger, UnionCharacter
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! ---------------------------------------------  FindAllCharacter  -----
  subroutine FindAllCharacter ( SET, IT, WHICH, HOW_MANY, RE_MAINDER, WHICH_NOT )
    ! Return which i of set[i] = it
    ! (case-sensitive, ignores trailing blanks, but alert to leading blanks)
    ! optionally, may return also:
    ! how many of them do
    ! which_not of them (which don't)
    ! and the re_mainder of the set != it
    ! e.g. given set = /(4, 3, 1, 2, 1, 3 )/ and it = 1
    ! produces which = /(3, 5)/, 
    !      which_not = /(1, 2, 4, 6)/, 
    !       how_many = 2,
    !     re_mainder = /(4, 3, 2, 3)/
    
    ! Note that which and which_not are arrays of array indices
    ! while re_mainder is an array of array elements
    
    ! This may be useful,
    ! e.g. in dump for reshaping an array to suppress any dims 
    ! that are identically 1
    
    ! Maybe this should be rewritten to use fraternal functions findfirst
    ! or findnext

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
  subroutine FindAllInteger ( SET, IT, WHICH, HOW_MANY, RE_MAINDER, WHICH_NOT )
    ! Return which i of set[i] = it
    ! optionally, may return also:
    ! how many of them do
    ! which_not of them (which don't)
    ! and the re_mainder of the set != it
    ! e.g. given set = /(4, 3, 1, 2, 1, 3 )/ and it = 1
    ! produces which = /(3, 5)/, 
    !      which_not = /(1, 2, 4, 6)/, 
    !       how_many = 2,
    !     re_mainder = /(4, 3, 2, 3)/
    
    ! Note that which and which_not are arrays of array indices
    ! while re_mainder is an array of array elements
    
    ! This may be useful,
    ! e.g. in dump for reshaping an array to suppress any dims 
    ! that are identically 1
    
    ! Maybe this should be rewritten to use fraternal functions findfirst
    ! or findnext

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
  subroutine FindAllLogical ( set, WHICH, HOW_MANY, WHICH_NOT )
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
  subroutine FindAllSubString ( SET, IT, WHICH, HOW_MANY, RE_MAINDER, WHICH_NOT )
    ! Return which substring i of set[i:i] = it
    ! Formal arguments
    character(len=*), intent(in)                 :: set
    character(len=1), intent(in)                 :: it
    integer, intent(out), dimension(:)           :: which
    integer, intent(out), optional               :: how_many
    character(len=*), intent(out), optional :: re_mainder
    integer, intent(out), dimension(:), optional :: which_not

    ! local variables
    integer :: i, i_which, i_re_mainder
    
    if ( len(set) < 1 .or. size(which) < 1 ) then
      if ( present(how_many) ) how_many = 0
      if ( present(re_mainder) ) re_mainder = ''
      return
    end if
    i_which = 0
    i_re_mainder = 0
    do i=1, len(set)
      if ( set(i:i) == it ) then
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

  end subroutine FindAllSubString

  ! -------------------------------------------  FindFirstCharacter  -----
  integer function FindFirstCharacter ( Set, Probe, Tol )
    ! Find the first element in the array Set that is equal to Probe
    ! (case-sensitive, ignores trailing blanks, but alert to leading blanks)
    character(len=*), dimension(:), intent(in) :: Set
    character(len=*), intent(in) :: Probe
    character(len=*), optional, intent(in) :: Tol ! Ignored; generic consistency

    ! Executable code
    do FindFirstCharacter = 1, size(set)
      if ( trim(set(FindFirstCharacter)) == trim(probe) ) return
    end do
    FindFirstCharacter = 0
  end function FindFirstCharacter

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

  ! ---------------------------------------------  FindIntersection  -----
  ! This family of routines finds the indices of the intersection between
  ! two similarly montonic sets
  ! e.g. given set1 = /(1, 3, 5, 7 )/, set2 = /(1, 2, 3, 4, 5)/
  ! produces which1 = /(1, 2, 3)/, 
  !      which2 = /(1, 3, 5)/, 
  ! and the optional arg      how_many = 3
  
  ! Note that which1 and which2 are arrays of array indices
  ! and that they will be monotonically increasing
  subroutine FindIntersectionInteger ( set1, set2, WHICH1, which2, &
    & HOW_MANY )
    ! Formal arguments
    integer, dimension(:), intent(in)     :: set1
    integer, dimension(:), intent(in)     :: set2
    integer, dimension(:), intent(out)    :: which1
    integer, dimension(:), intent(out)    :: which2
    integer, optional, intent(out)        :: how_many
    ! Internal variables
    integer                               :: myTol
    integer :: i1, i2
    integer :: my_how_many
    ! Executable
    myTol = 0
    include 'FindIntersection.f9h'
  end subroutine FindIntersectionInteger

  subroutine FindIntersectionReal ( set1, set2, WHICH1, which2, &
    & HOW_MANY, tol )
    ! Formal arguments
    real, dimension(:), intent(in)        :: set1
    real, dimension(:), intent(in)        :: set2
    integer, dimension(:), intent(out)    :: which1
    integer, dimension(:), intent(out)    :: which2
    integer, optional, intent(out)        :: how_many
    real, intent(in), optional            :: Tol
    ! Internal variables
    real                                  :: myTol
    integer :: i1, i2
    integer :: my_how_many
    ! Executable
    myTol = 0.
    if ( present(tol) ) myTol = tol
    include 'FindIntersection.f9h'
  end subroutine FindIntersectionReal

  subroutine FindIntersectionDouble ( set1, set2, WHICH1, which2, &
    & HOW_MANY, tol )
    ! Formal arguments
    double precision, dimension(:), intent(in)   :: set1
    double precision, dimension(:), intent(in)   :: set2
    integer, dimension(:), intent(out)           :: which1
    integer, dimension(:), intent(out)           :: which2
    integer, optional, intent(out)               :: how_many
    double precision, intent(in), optional       :: Tol
    ! Internal variables
    double precision                             :: myTol
    integer :: i1, i2
    integer :: my_how_many
    ! Executable
    myTol = 0.d0
    if ( present(tol) ) myTol = tol
    include 'FindIntersection.f9h'
  end subroutine FindIntersectionDouble

  subroutine FindIntersectionCharacter ( set1, set2, WHICH1, which2, &
    & HOW_MANY )
    ! Formal arguments
    character(len=*), dimension(:), intent(in)  :: set1
    character(len=*), dimension(:), intent(in)  :: set2
    integer, dimension(:), intent(out)          :: which1
    integer, dimension(:), intent(out)          :: which2
    integer, optional, intent(out)              :: how_many
    ! Internal variables
    character(len=1)                               :: myTol
    integer :: i1, i2
    integer :: my_how_many
    ! Executable
    myTol = ''
    include 'FindIntersection.f9h'
  end subroutine FindIntersectionCharacter

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

  ! --------------------------------------------  FindNextCharacter  -----
  integer function FindNextCharacter ( Set, Probe, Current, Wrap, Repeat )
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

    ! Local variables
    integer :: I                        ! Loop counter
    logical :: myWrap
    logical :: myRepeat

    ! Executable code
    myWrap = .false.
    if ( present(wrap) ) myWrap = wrap
    myRepeat = .false.
    if ( present(repeat) ) myRepeat = repeat
    FindNextCharacter = 0
    ! We'll assume you gave us valid args; otherwise return 0
    if ( current < 1 .or. current > size(set)) return
    if ( set(current) /= probe ) return
    ! Now check for current already at end of array
    if ( current < size(set) ) then
      do i = current+1, size(set)
        if ( trim(set(i)) == trim(probe) ) then
          FindNextCharacter = i
          return
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
    if ( set(current) /= probe ) return
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
    if ( .not. condition(current)) return
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
    if ( (set(current:current) == probe) .eqv. myReverse ) return
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

  ! --------------------------------------------------  Intersect  -----
  ! This family of functions returns TRUE if the two sets intersect,
  ! i.e. have at least one element in common
  logical function IntersectInteger ( A, B )
  ! Return true if the integer arrays A and B have a common element

    integer, intent(in) :: A(:), B(:)

    integer :: I, J

    ! Put the longer loop as the inner one
    if ( size(a) < size(b) ) then
      do i = 1, size(a)
        do j = 1, size(b)
          if ( a(i) == b(j) ) then
            IntersectInteger = .true.
            return
          end if
        end do
      end do
    else
      do i = 1, size(b)
        do j = 1, size(a)
          if ( a(j) == b(i) ) then
            IntersectInteger = .true.
            return
          end if
        end do
      end do
    end if
    IntersectInteger = .false.

  end function IntersectInteger

  logical function IntersectCharacter ( A, B )
  ! Return true if the character arrays A and B have a common element

    character(len=*), intent(in) :: A(:), B(:)

    integer :: I, J

    ! Put the longer loop as the inner one
    if ( size(a) < size(b) ) then
      do i = 1, size(a)
        do j = 1, size(b)
          if ( a(i) == b(j) ) then
            IntersectCharacter = .true.
            return
          end if
        end do
      end do
    else
      do i = 1, size(b)
        do j = 1, size(a)
          if ( a(j) == b(i) ) then
            IntersectCharacter = .true.
            return
          end if
        end do
      end do
    end if
    IntersectCharacter = .false.

  end function IntersectCharacter

  ! -----------------------------------------------  Intersection  -----
  ! Compute the intersection C of the sets A and B, each represented by
  ! arrays of integers, characters, or reals
  
  ! options, if present, can modify this behavior: 
  !   char    effect
  !  ------   ------
  !    r      The set of elements in a or b but not in both
  !            which is Union(A, B) - Intersection(A, B)
  !    c      The set of elements in a but not in b
  !            which is A - B, or as our notation has it, RC ( B, A )
  ! The 'c' option exists only to allow the family of functions
  ! to be reused internally for finding the relative complement of two sets
  
  function IntersectionInteger ( A, B, options ) result ( C )
    ! A faster algorithm is used if we're not reversing

    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_ERROR
    use Sort_M, only: SORT

    integer, intent(in) :: A(:), B(:)
    integer, pointer :: C(:) ! Intent(out) -- nullified and then allocated here
    ! logical, optional, intent(in) :: reverse
    character(len=*), optional, intent(in) :: options

    integer :: size_c, status
    integer :: I, J, K, Stat, TA(size(a)), TB(size(b)), TC(size(a)+size(b))
    logical :: myComplement
    logical :: myReverse
    logical :: stdIntersection

    ! Executable
    myComplement = .false.
    myReverse = .false.
    if ( present(options) ) myComplement = index( options, 'c' ) > 0
    if ( present(options) ) myReverse = index( options, 'r' ) > 0
    stdIntersection = .not. (myComplement .or. myReverse )
    if ( .not. stdIntersection ) then
    include 'Intersection.f9h'
    else
      ta = a
      tb = b
      call sort ( ta, 1, size(ta) )
      call sort ( tb, 1, size(tb) )

      i = 1; j=1; k=0
      do while ( i <= size(ta) .and. j <= size(tb) )
        if ( ta(i) == tb(j) ) then
          tc(k+1) = ta(i)
          i = i + 1; j = j + 1; k = k + 1
        else if ( ta(i) < tb(j) ) then
          i = i + 1
        else
          j = j + 1
        end if
      end do

      nullify ( c )
      allocate ( c(k), stat=stat )
      if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        MLSMSG_Allocate // 'C in IntersectionInteger' )
      c = tc(:k)
    end if
  end function IntersectionInteger

  function IntersectionCharacter ( A, B, options ) result ( C )
    ! method:
    ! Go though a, checking for each element whether a match is found in (b)
    ! If  so found, add the element
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_ERROR

    character(len=*), dimension(:), intent(in) :: A(:), B(:)
    character(len=len(a)), dimension(:), pointer :: C(:) ! Intent(out) -- nullified and then allocated here
    ! logical, optional, intent(in) :: reverse
    character(len=*), optional, intent(in) :: options
    ! Local variables
    integer :: i, j, size_c, status
    character(len=len(a)), dimension(size(a)+size(b)) :: TC
    logical :: myComplement
    logical :: myReverse
    logical :: stdIntersection
    
    include 'Intersection.f9h'
  end function IntersectionCharacter

  function IntersectionDouble ( A, B, options ) result ( C )
    ! method:
    ! Go though a, checking for each element whether a match is found in (b)
    ! If  so found, add the element
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_ERROR

    double precision, dimension(:), intent(in) :: A(:), B(:)
    double precision, dimension(:), pointer :: C(:) ! Intent(out) -- nullified and then allocated here
    ! logical, optional, intent(in) :: reverse
    character(len=*), optional, intent(in) :: options
    ! Local variables
    integer :: i, j, size_c, status
    double precision, dimension(size(a)+size(b)) :: TC
    logical :: myComplement
    logical :: myReverse
    logical :: stdIntersection
    
    include 'Intersection.f9h'
  end function IntersectionDouble

  function IntersectionReal ( A, B, options ) result ( C )
    ! method:
    ! Go though a, checking for each element whether a match is found in (b)
    ! If  so found, add the element
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_ERROR

    real, dimension(:), intent(in) :: A(:), B(:)
    real, dimension(:), pointer :: C(:) ! Intent(out) -- nullified and then allocated here
    ! logical, optional, intent(in) :: reverse
    character(len=*), optional, intent(in) :: options
    ! Local variables
    integer :: i, j, size_c, status
    real, dimension(size(a)+size(b)) :: TC
    logical :: myComplement
    logical :: myReverse
    logical :: stdIntersection
    
    include 'Intersection.f9h'
  end function IntersectionReal

  ! --------------------------------------------------  IsProperSet  -----
  ! This family of functions returns TRUE if A is a proper set; i.e.
  ! each element of A is unique
  logical function IsProperSetInteger ( A )

    integer, intent(in) :: A(:)

    integer :: nUnique
    integer, dimension(size(A)) :: B

    call FindUnique( A, B, NUnique )
    IsProperSetInteger = (NUnique == size(A))

  end function IsProperSetInteger

  logical function IsProperSetCharacter ( A )
  ! Return true if the character arrays A and B have a common element

    character(len=*), intent(in) :: A(:)

    integer :: nUnique
    character(len=len(A(1))), dimension(size(A)) :: B

    call FindUnique( A, B, NUnique )
    IsProperSetCharacter = (NUnique == size(A))

  end function IsProperSetCharacter

  ! --------------------------------------------------  IsProperSubset  -----
  ! This family of functions returns TRUE if the A is a proper subset of B
  ! element of A is in B, but there is at least one element of B not in A
  logical function IsProperSubsetInteger ( A, B )

    integer, intent(in) :: A(:), B(:)

    integer :: I, J

    ! Is B the longer array? If not, we're already done
    IsProperSubsetInteger = .false.
    if ( size(a) < size(b) ) IsProperSubsetInteger = IsSubsetInteger ( A, B )

  end function IsProperSubsetInteger

  logical function IsProperSubsetCharacter ( A, B )
  ! Return true if the character arrays A and B have a common element

    character(len=*), intent(in) :: A(:), B(:)

    integer :: I, J

    ! Is B the longer array? If not, we're already done
    IsProperSubsetCharacter = .false.
    if ( size(a) < size(b) ) IsProperSubsetCharacter = IsSubsetCharacter ( A, B )
  end function IsProperSubsetCharacter

  ! --------------------------------------------------  IsSubset  -----
  ! This family of functions returns TRUE if the A is a subset of B
  ! i.e., each element of A is in B
  logical function IsSubsetInteger ( A, B )

    integer, intent(in) :: A(:), B(:)

    integer :: I, J

    IsSubsetInteger = .false.
    do i = 1, size(a)
      j = FindFirst (b, a(i))
      if ( j < 1 ) return
    end do
    IsSubsetInteger = .true.

  end function IsSubsetInteger

  logical function IsSubsetCharacter ( A, B )
  ! Return true if the character arrays A and B have a common element

    character(len=*), intent(in) :: A(:), B(:)

    integer :: I, J

    IsSubsetCharacter = .false.
    do i = 1, size(a)
      j = FindFirst (b, a(i))
      if ( j < 1 ) return
    end do
    IsSubsetCharacter = .true.
  end function IsSubsetCharacter

  ! -----------------------------------------------  RelativeComplement  -----
  ! Compute the RelativeComplement C of the sets A and B, each represented by
  ! arrays of integers or characters
  !
  ! Example:
  ! Let A = {a, b, c} and B = {a, b, d, e, f, g}
  ! Then C = {d, e, f, g} which are all the elements of B except those in A
  ! If B is the "Universal Set" then this would be the complement of A
  ! in which case there should not be any elements in A not also found in B

  function RelativeComplementInteger ( A, B ) result ( C )

    integer, intent(in) :: A(:), B(:)
    integer, pointer :: C(:) ! Intent(out) -- nullified and then allocated here

    C => Intersection ( B, A, options = 'c' )
  end function RelativeComplementInteger

  function RelativeComplementCharacter ( A, B ) result ( C )
    ! method:
    ! Go though b, checking for each element whether a match is found in (a)
    ! If  not, add the element

    character(len=*), dimension(:), intent(in) :: A(:), B(:)
    character(len=len(a)), dimension(:), pointer :: C(:) ! Intent(out) -- nullified and then allocated here

    C => Intersection ( B, A, options = 'c' )
  end function RelativeComplementCharacter

  function RelativeComplementDouble ( A, B ) result ( C )

    double precision, intent(in) :: A(:), B(:)
    double precision, pointer :: C(:) ! Intent(out) -- nullified and then allocated here

    C => Intersection ( B, A, options = 'c' )
  end function RelativeComplementDouble

  function RelativeComplementReal ( A, B ) result ( C )

    real, intent(in) :: A(:), B(:)
    real, pointer :: C(:) ! Intent(out) -- nullified and then allocated here

    C => Intersection ( B, A, options = 'c' )
  end function RelativeComplementReal

  ! ------------------------------------------------------  Union  -----
  ! Compute the union C of the sets A and B, each represented by
  ! arrays of integers, characters
  function UnionInteger ( A, B ) result ( C )

    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_ERROR
    use Sort_M, only: SORT

    integer, intent(in) :: A(:), B(:)
    integer, pointer :: C(:) ! Intent(out) -- nullified and then allocated here

    integer :: I, J, Stat, T(size(a)+size(b))

    t(1:size(a)) = a
    t(size(a)+1:size(t)) = b
    call sort ( t, 1, size(t) )

    ! remove duplicates
    i = 0; j = 0
    do while ( i < size(t) .and. j < size(t) )
      i = i + 1; j = j + 1
      t(i) = t(j)
      do while ( j < size(t) )
        if ( t(j+1) /= t(i) ) exit
        j = j + 1
      end do
    end do

    allocate ( c(i), stat=stat )
    if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      MLSMSG_Allocate // 'C in Intersection' )
    c = t(:i)

  end function UnionInteger

  function UnionCharacter ( A, B ) result ( C )
    ! Method:
    ! 1st go through a, adding all non-repeated elements
    ! Then go through b, adding any that haven't been added already
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_ERROR

    character(len=*), dimension(:), intent(in) :: A(:), B(:)
    character(len=len(a)), dimension(:), pointer :: C(:) ! Intent(out) -- nullified and then allocated here
    ! Local variables
    integer :: i, j, size_c, status
    character(len=len(a)), dimension(size(a)+size(b)) :: TC
    
    ! Executable
    ! print *, 'size(a): ', size(a)    
    ! print *, 'a: ', a(1:size(a))
    size_c = 0
    do i=1, size(a)
      ! Don't redo a repeated element
      if ( i > 1 ) then
        j = findFirst( a(:i-1), a(i) )
        if ( j > 0 ) cycle
      end if
      size_c = size_c + 1
      TC(size_c) = a(i)
    end do

    ! print *, 'size(b): ', size(b)    
    ! print *, 'b: ', b(1:size(b))
    do i=1, size(b)
      ! Don't redo an already-added element
      if ( size_c > 0 ) then
        j = findFirst( tc(:size_c), b(i) )
        if ( j > 0 ) cycle
      end if
      size_c = size_c + 1
      TC(size_c) = b(i)
    end do
    ! print *, 'size(c): ', size_c    
    ! print *, 'tc: ', tc(1:size_c)
    nullify(c)
    if ( size_c < 1 ) return
    allocate ( c(size_c), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      MLSMSG_Allocate // 'C in UnionCharacter' )
    c = tc(:size_c)
  end function UnionCharacter

  ! --------------------------------------------------  UnionSize  -----
  integer function UnionSize ( A, B )
  ! Compute the size of the union of the sets A and B, each represented by
  ! arrays of integers.

    use Sort_M, only: SORT

    integer, intent(in) :: A(:), B(:)

    integer :: I, J, T(size(a)+size(b))

    t(1:size(a)) = a
    t(size(a)+1:size(t)) = b
    call sort ( t, 1, size(t) )

    ! remove duplicates
    i = 0; j = 0
    do while ( i < size(t) .and. j < size(t) )
      i = i + 1; j = j + 1
      t(i) = t(j)
      do while ( j < size(t) )
        if ( t(j+1) /= t(i) ) exit
        j = j + 1
      end do
    end do

    unionSize = i

  end function UnionSize

! =====     Private Procedures     =====================================

  subroutine FindLongestStretch( which, how_many, range )
    ! Given an array of integers, return the longest
    ! stretch of consecutive ones
    integer, dimension(:), intent(in)  :: which
    integer, intent(in)                :: how_many
    integer, dimension(2), intent(out) :: range
    ! Internal variables
    integer :: firstWhich
    integer :: i
    logical :: inSequence
    integer :: lastWhich
    integer :: stretch
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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSSets

! $Log$
! Revision 2.25  2012/12/04 00:10:57  pwagner
! Changed api to Intersection; improced module toc and api
!
! Revision 2.24  2012/08/07 17:59:47  pwagner
! FindLast.. now takes optional arg Reverse
!
! Revision 2.23  2011/04/28 22:43:29  vsnyder
! Regularize declarations in FindFirst... routines
!
! Revision 2.22  2011/02/18 18:03:14  pwagner
! Added functions to check set A is proper or that A is a [proper] subset of B
!
! Revision 2.21  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.20  2008/11/24 19:39:19  pwagner
! Added RelativeComplement function
!
! Revision 2.19  2008/06/18 20:45:25  pwagner
! FindUnique can now take real, double Sets
!
! Revision 2.18  2008/02/22 21:25:53  pwagner
! FindFirst can handle periods
!
! Revision 2.17  2007/12/19 01:27:42  pwagner
! Added FindUnique to reduce input Set to its unique members
!
! Revision 2.16  2007/10/09 00:30:50  pwagner
! Added FindLongestRange procedures
!
! Revision 2.15  2007/09/20 18:38:31  pwagner
! Added substring versions of of Find(-First, Last, Next)
!
! Revision 2.14  2007/02/26 23:54:39  pwagner
! Added FindIntersection; FindFirst can search through arrays of any numerical type
!
! Revision 2.13  2006/08/12 00:07:43  pwagner
! Added 2d findFirst, Last
!
! Revision 2.12  2006/07/24 20:36:22  pwagner
! Union, Inersection may take character arrays
!
! Revision 2.11  2006/01/14 00:51:39  pwagner
! Added FindLast functions
!
! Revision 2.10  2005/06/07 00:49:24  vsnyder
! Status quo ante 2.9 -- FindFirst means _first_ not _last_
!
! Revision 2.9  2005/06/04 00:32:28  vsnyder
! New copyright, move Id to not_used_here, simplify FindFirst...
!
! Revision 2.8  2004/07/02 01:34:11  vsnyder
! Get rid of unused stuff
!
! Revision 2.7  2004/06/16 22:14:16  pwagner
! character arrays can bow be args of find... methods
!
! Revision 2.6  2004/06/16 01:24:38  vsnyder
! Add UnionSize
!
! Revision 2.5  2004/06/11 20:03:01  vsnyder
! Cannonball polishing
!
! Revision 2.4  2004/06/11 20:02:14  vsnyder
! Add Intersect function
!
! Revision 2.3  2004/06/11 19:03:35  pwagner
! Added FindAll
!
! Revision 2.2  2004/06/10 20:03:14  vsnyder
! Don't use Allocate_Test -- results in USE cycle
!
! Revision 2.1  2004/06/10 00:12:22  vsnyder
! Initial commit
!
