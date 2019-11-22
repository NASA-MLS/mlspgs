! Copyright 2013, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MLSSets

! Various abstract operations on sets and sequences

! Note:
! A set in mathematics consists of elements
! Each element in a set appears exactly once, i.e. each is unique

! We represent sets with arrays
! We usually do not check whether the elements are in fact unique
! (you can use IsProperSet to check or FindUnique to ensure this)
! so some operations may not work properly if you supply arrays with
! duplicate elements.

! You may also call Union with any args, proper or not,
! confident that it will return only unique elements.

! A point about null sets:
! Not all operations have been checked thoroughly for correctness
! when one or more arguments are null sets, i.e. 0-sized
! IsProperSet, IsProperSubset, and IsSubset do give correct results:
! The null set is
!  (1) a proper set,
!  (2) a proper subset of any set except itself
!  (3) a subset of itself, but not a proper subset of itself

! The order in which elements appear is immaterial for sets.
! For sequences, however, the order is essential
! Many procedures in this module treat sets with little regard for order
! FindFirst, -Next, -Last, and -Longest functions must
! pay attention to order, however, so properly speaking
! they treat sequences, not sets

! Furthermore, a mathematical set's elements need not all be of the same type
! while a sequence's elements are necessarily homogeneous. From this
! standpoint, we might say that all the procedures in this module
! treat sequences, not sets, even when order is immaterial
  use MLSFinds, only: FindFirst, FindUnique

  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (subroutines and functions)
! AreSetsEqual  Check that two sets are equal, i.e. have the same elements
! FindIntersection
!               Compute indices of elements in intersection between two sets
! Intersect     Return true if two sets share at least one common element
! Intersection  Compute intersection between two sets
! IsProperSet   Check that each element is unique
! IsProperSubset
!               Check that each element in A is within larger B
! IsSubset      Check that each element in A is within B
! RelativeComplement
!               Compute complement of set a relative to set b
!               i.e., all elements in b except those in a
! Union         Compute the union between two sets
! UnionSize     Compute the size that the union between two sets would have.
! === (end of toc) ===

! === (start of api) ===
! log AreSetsEqual ( set A, set B )
! FindIntersection ( set1(:), set2(:), int which1(:), int which2(:),
!      [int how_many] )
! log Intersect ( set A, set B )
! set Intersection ( set A, set B, [char* options] )
! log isProperSet ( set A )
! log isProperSubset ( set A, set B )
! log isSubset ( set A, set B )
! set RelativeComplement ( set A, set B )
! set Union ( set A, set B, [int how_many] )
! int UnionSize ( set A, set B )

! Notes:
! A set is represented either
! as an array of integers, floats, or an array of characters; e.g.
!   logical Intersect ( set a(:), set b(:) )
!   set *Intersection ( set a(:), set b(:) )
!   logical IsProperSet ( set a(:) )
!   logical Is[Proper]Subset ( set a(:), set b(:) )
!   set *RelativeComplement ( set a(:), set b(:) )
!   set *Union ( set a(:), set b(:) )
!   int UnionSize ( set a(:), set b(:) )
! === (end of api) ===

  public :: AreSetsEqual, FindIntersection, &
    & Intersect, Intersection, IsProperSet, IsProperSubset, IsSubset, &
    & RelativeComplement, &
    & Union, Unionsize

  interface AreSetsEqual
    module procedure AreSetsEqualInteger, AreSetsEqualCharacter
  end interface

  interface FindIntersection
    module procedure FindIntersectionInteger, FindIntersectionReal, &
      & FindIntersectionDouble, FindIntersectionCharacter
  end interface

  interface Intersect
    module procedure IntersectInteger, IntersectCharacter
  end interface

  interface Intersection
    module procedure IntersectionInteger, IntersectionCharacter
    module procedure IntersectionReal, IntersectionDouble, &
      & IntersectionVector
  end interface

  interface IsProperSet
    module procedure IsProperSetInteger, IsProperSetCharacter
  end interface

  interface IsProperSubset
    module procedure IsProperSubsetInteger, IsProperSubsetCharacter
  end interface

  interface IsSubset
    module procedure IsSubsetInteger, IsSubsetCharacter
  end interface

  interface RelativeComplement
    module procedure RelativeComplementInteger, RelativeComplementCharacter
    module procedure RelativeComplementReal, RelativeComplementDouble
  end interface

  interface Union
    module procedure UnionInteger, UnionCharacter, UnionVector
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! --------------------------------------------------  AreSetsEqual  -----
  ! This family of functions returns TRUE if the sets A and B are the same
  ! i.e., each element of A is in B and vice versa
  logical function AreSetsEqualInteger ( A, B )
    integer, intent(in) :: A(:), B(:)
    AreSetsEqualInteger = IsSubset( A, B ) .and. IsSubset( B, A )
  end function AreSetsEqualInteger

  logical function AreSetsEqualCharacter ( A, B )
    character(len=*), intent(in) :: A(:), B(:)
    AreSetsEqualCharacter = IsSubset( A, B ) .and. IsSubset( B, A )
  end function AreSetsEqualCharacter

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
  ! Compute the intersection C between the sets A and B, each represented by
  ! arrays of integers, characters, or reals

  ! options, if present, can modify this behavior:
  !   char    effect
  !  ------   ------
  !    r      The set of elements in a or b but not in both
  !            which is Union(A, B) - Intersection(A, B)
  !    c      The set of elements in a but not in b
  !            which is A - B, or as our notation has it, RC ( B, A )
  ! The 'c' option exists only to allow the family of functions
  ! to be reused internally for finding the relative complement between two sets

  function IntersectionInteger ( A, B, options ) result ( C )
    ! A faster algorithm is used if we're not reversing

    use Sort_M, only: Sort

    integer, dimension(:), intent(in)      :: A, B
    integer, dimension(:), allocatable     :: C ! Intent(out) -- allocated here
    character(len=*), optional, intent(in) :: options

    integer :: size_c
    integer :: I, J, K, TA(size(a)), TB(size(b)), TC(size(a)+size(b))
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

      c = tc(:k)
    end if
  end function IntersectionInteger

  function IntersectionCharacter ( A, B, options ) result ( C )
    ! method:
    ! Go though a, checking for each element whether a match is found in (b)
    ! If  so found, add the element

    character(len=*), dimension(:), intent(in)   :: A, B
    character(len=len(a)), allocatable :: C(:) ! Intent(out) -- allocated here
    character(len=*), optional, intent(in)       :: options
    ! Local variables
    integer :: i, j, size_c
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

    double precision, dimension(:), intent(in)  :: A, B
    double precision, dimension(:), allocatable :: C ! Intent(out) -- allocated here
    character(len=*), optional, intent(in)      :: options
    ! Local variables
    integer :: i, j, size_c
    double precision, dimension(size(a)+size(b)) :: TC
    logical :: myComplement
    logical :: myReverse
    logical :: stdIntersection

    include 'Intersection.f9h'
  end function IntersectionDouble

  function IntersectionReal ( A, B, options ) result ( C )
    ! method:
    ! Go though A, checking for each element whether a match is found in (B)
    ! If  so found, add the element

    real, dimension(:), intent(in)         :: A, B
    real, dimension(:), allocatable        :: C ! Intent(out) -- allocated here
    character(len=*), optional, intent(in) :: options
    ! Local variables
    integer :: i, j, size_c
    real, dimension(size(a)+size(b)) :: TC
    logical :: myComplement
    logical :: myReverse
    logical :: stdIntersection

    include 'Intersection.f9h'
  end function IntersectionReal

  function IntersectionVector ( A, B, options ) result ( C )
    ! A vector is represented as an array of integers (1d)
    ! Thus an array of vectors is a 2d array of integers
    ! A and B must have the same lower and upper bounds for their first index

    integer, dimension(:,:), intent(in) :: A, B
    integer, allocatable :: C(:,:) ! Intent(out) -- allocated here
    character(len=*), optional, intent(in) :: options  ! then allocated here
    ! Local variables
    integer :: i, j, size_c
    integer, dimension(size(a,1),size(a,2)+size(b,2)) :: TC
    logical :: myComplement
    logical :: myReverse
    logical :: stdIntersection

    ! Executable
    myComplement = .false.
    myReverse = .false.
    if ( present(options) ) myComplement = index( options, 'c' ) > 0
    if ( present(options) ) myReverse = index( options, 'r' ) > 0
    stdIntersection = .not. (myComplement .or. myReverse )
    size_c = 0
    do i=1, size(a,2)
      ! Don't redo a repeated element
      if ( i > 1 ) then
        j = findFirst( a(:,:i-1), a(:,i) )
        if ( j > 0 ) cycle
      endif
      j = findFirst( b, a(:,i) )
      if ( (stdIntersection .and. j > 0) .or. &
        & (.not. stdIntersection .and. j < 1) ) then
        size_c = size_c + 1
        TC(:,size_c) = a(:,i)
      endif
    enddo
    ! Now switch to looping over B
    if ( myReverse ) then
      do i=1, size(b,2)
        ! Don't redo a repeated element
        if ( i > 1 ) then
          j = findFirst( b(:,:i-1), b(:,i) )
          if ( j > 0 ) cycle
        endif
        j = findFirst( a, b(:,i) )
        if ( j < 1 ) then
          size_c = size_c + 1
          TC(:,size_c) = b(:,i)
        endif
      enddo
    endif
    ! print *, 'size(c): ', size_c
    ! print *, 'tc: ', tc(1:size_c)
    if ( size_c < 1 ) return
    c = tc(:,:size_c)
  end function IntersectionVector

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

    ! Is B the longer array? If not, we're already done
    IsProperSubsetInteger = .false.
    if ( size(a) < size(b) ) IsProperSubsetInteger = IsSubsetInteger ( A, B )

  end function IsProperSubsetInteger

  logical function IsProperSubsetCharacter ( A, B )
  ! Return true if the character arrays A and B have a common element

    character(len=*), intent(in) :: A(:), B(:)

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
    integer, allocatable :: C(:) ! Intent(out) -- allocated here

!   call move_alloc ( Intersection ( B, A, options = 'c' ), c ) doesn't work
    C = Intersection ( B, A, options = 'c' )
  end function RelativeComplementInteger

  function RelativeComplementCharacter ( A, B ) result ( C )
    ! method:
    ! Go though b, checking for each element whether a match is found in (a)
    ! If  not, add the element

    character(len=*), dimension(:), intent(in) :: A, B
    character(len=len(a)), allocatable :: C(:) ! Intent(out) -- allocated here

!   call move_alloc ( Intersection ( B, A, options = 'c' ), c ) doesn't work
    C = Intersection ( B, A, options = 'c' )
  end function RelativeComplementCharacter

  function RelativeComplementDouble ( A, B ) result ( C )

    double precision, intent(in) :: A(:), B(:)
    double precision, allocatable :: C(:) ! Intent(out) - allocated here

!   call move_alloc ( Intersection ( B, A, options = 'c' ), c ) doesn't work
    C = Intersection ( B, A, options = 'c' )
  end function RelativeComplementDouble

  function RelativeComplementReal ( A, B ) result ( C )

    real, intent(in) :: A(:), B(:)
    real, allocatable :: C(:) ! Intent(out) -- allocated here

!   call move_alloc ( Intersection ( B, A, options = 'c' ), c ) doesn't work
    C = Intersection ( B, A, options = 'c' )
  end function RelativeComplementReal

  ! ------------------------------------------------------  Union  -----
  ! Compute the union C of the sets A and B, each represented by
  ! arrays of integers, characters
  function UnionInteger ( A, B, how_many ) result ( C )

    use Sort_M, only: Sort

    integer, intent(in)                  :: A(:), B(:)
    integer, allocatable                 :: C(:) ! Intent(out) -- allocated here
    integer, optional, intent(out)       :: how_many

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

    c = t(:i)
    if ( present(how_many) ) how_many = i

  end function UnionInteger

  function UnionCharacter ( A, B, how_many ) result ( C )
    ! Method:
    ! 1st go through a, adding all non-repeated elements
    ! Then go through b, adding any that haven't been added already

    character(len=*), dimension(:), intent(in) :: A, B
    character(len=len(a)), allocatable         :: C(:) ! Intent(out) -- allocated here
    integer, optional, intent(out)             :: how_many
    ! Local variables
    integer :: i, j, size_c
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
    c = tc(:size_c)
    if ( present(how_many) ) how_many = size_c
  end function UnionCharacter

  ! --------------------------------------------------  UnionSize  -----
  integer function UnionSize ( A, B )
  ! Compute the size of the union of the sets A and B, each represented by
  ! arrays of integers.

    use Sort_M, only: Sort

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

  function UnionVector ( A, B ) result ( C )
    ! A vector is represented as an array of integers (1d)
    ! Thus an array of vectors is a 2d array of integers
    ! A and B must have the same lower and upper bounds for their first index

    integer, dimension(:,:), intent(in) :: A, B
    integer, allocatable :: C(:,:) ! Intent(out) -- allocated here
    ! Local variables
    integer :: i, j, size_c
    integer, dimension(size(a,1),size(a,2)+size(b,2)) :: TC

    ! Executable
    ! print *, 'size(a): ', size(a)
    ! print *, 'a: ', a(1:size(a))
    size_c = 0
    do i=1, size(a,2)
      ! Don't redo a repeated element
      if ( i > 1 ) then
        j = findFirst( a(:,:i-1), a(:,i) )
        if ( j > 0 ) cycle
      end if
      size_c = size_c + 1
      TC(:,size_c) = a(:,i)
    end do

    ! print *, 'size(b): ', size(b)
    ! print *, 'b: ', b(1:size(b))
    do i=1, size(b,2)
      ! Don't redo an already-added element
      if ( size_c > 0 ) then
        j = findFirst( tc(:,:size_c), b(:,i) )
        if ( j > 0 ) cycle
      end if
      size_c = size_c + 1
      TC(:,size_c) = b(:,i)
    end do
    ! print *, 'size(c): ', size_c
    ! print *, 'tc: ', tc(1:size_c)
    c = tc(:,:size_c)
  end function UnionVector

! =====     Private Procedures     =====================================

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
! Revision 2.40  2019/11/22 00:35:11  pwagner
! Some of the procedures may opttionally return the size of the result
!
! Revision 2.39  2018/08/02 23:09:10  pwagner
! Added AreSetsEqual
!
! Revision 2.38  2018/04/18 23:15:44  vsnyder
! Remove declarations of unused variables
!
! Revision 2.37  2018/04/18 22:38:07  vsnyder
! Make results of IntersectionDouble and IntersectionReal allocatable instead
! of pointer, so they don't leak.  Return zero-size instead of unallocated
! array for zero-size set result.
!
! Revision 2.36  2017/11/02 00:07:52  pwagner
! Intersection and Union can now take vectors as args
!
! Revision 2.35  2017/03/10 01:00:02  vsnyder
! Remove declarations of unused variables
!
! Revision 2.34  2017/03/10 00:59:01  vsnyder
! Return an empty set instead of undefined result from UnionCharacter
!
! Revision 2.33  2017/03/10 00:38:10  vsnyder
! Make results of Intersection, RelativeComplement, Union allocatable
!
! Revision 2.32  2014/08/06 23:19:11  vsnyder
! Remove USE for unreferenced identifier FINDALL
!
! Revision 2.31  2013/08/12 23:46:20  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.30  2013/07/30 23:25:31  pwagner
! FindFirst now treats integer vectors and runs, too
!
! Revision 2.29  2013/06/18 17:58:29  pwagner
! Removed unused stuff; corrected syntax NAG found Questionable
!
! Revision 2.28  2013/06/17 21:37:15  pwagner
! FindNext no longer bails when current is not a solution
!
! Revision 2.27  2013/01/16 22:17:20  pwagner
! Fixed a syntax error only NAG noticed
!
! Revision 2.26  2013/01/15 18:54:44  pwagner
! Added FindLast for reals, doubles
!
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
