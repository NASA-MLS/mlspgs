! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MLSSets

! Various operations on sets

  implicit NONE
  private
  public :: FindAll, FindFirst, FindNext, Intersect, Intersection, Union, &
    & UnionSize

  interface FindFirst
    module procedure FindFirstInteger, FindFirstLogical, FindFirstCharacter
  end interface

  interface FindNext
    module procedure FindNextInteger, FindNextLogical, FindNextCharacter
  end interface

  interface FindAll
    module procedure FindAllInteger, FindAllLogical, FindAllCharacter
  end interface

! === (start of toc) ===                                                 
!     c o n t e n t s                                                    
!     - - - - - - - -                                                    

!     (subroutines and functions)
! FindAll       Find all logicals in the array that are true, or all the
!               integers in the array equal to the probe
! FindFirst     Find the first logical in the array that is true, or the
!               first integer in the array equal to the probe
! FindNext      Find the next logical in the array that is true, or the
!               next integer in the array equal to the probe
! Intersect     Return true if two sets represented by arrays of integers have
!               a common element
! Intersection  Compute intersection of two sets, represented as arrays of integers
! Union         Compute union of two sets, represented as arrays of integers
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
! int FindFirstInteger (int set(:), int probe)      
! int FindFirstLogical (log condition(:))      
! int FindNextCharacter (char* set(:), char* probe, int current, {log wrap}, {log repeat})      
! int FindNextInteger (int set(:), int probe, int current, {log wrap}, {log repeat})      
! int FindNextLogical (log condition(:), int current, {log wrap}, {log repeat})
! logical Intersect ( int a(:), int b(:) )
! int *Intersection ( int a(:), int b(:) )
! int *Union ( int a(:), int b(:) )
! int UnionSize ( int a(:), int b(:) )
! === (end of api) ===

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
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

  ! -------------------------------------------  FindFirstCharacter  -----
  integer function FindFirstCharacter ( Set, Probe )
    ! Find the first element in the array Set that is equal to Probe
    ! (case-sensitive, ignores trailing blanks, but alert to leading blanks)
    character(len=*), dimension(:), intent(in) :: Set
    character(len=*), intent(in) :: Probe

    ! Executable code
    do FindFirstCharacter = 1, size(set)
      if ( trim(set(FindFirstCharacter)) == trim(probe) ) return
    end do
    FindFirstCharacter = 0
  end function FindFirstCharacter

  ! -------------------------------------------  FindFirstInteger  -----
  integer function FindFirstInteger ( Set, Probe )
    ! Find the first element in the array Set that is equal to Probe
    integer, dimension(:), intent(in) :: Set
    integer, intent(in) :: Probe

    ! Executable code
    do FindFirstInteger = 1, size(set)
      if ( set(FindFirstInteger) == probe ) return
    end do
    FindFirstInteger = 0
  end function FindFirstInteger

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

  ! --------------------------------------------------  Intersect  -----
  logical function Intersect ( A, B )
  ! Return true if the integer arrays A and B have a common element

    integer, intent(in) :: A(:), B(:)

    integer :: I, J

    ! Put the longer loop as the inner one
    if ( size(a) < size(b) ) then
      do i = 1, size(a)
        do j = 1, size(b)
          if ( a(i) == b(j) ) then
            intersect = .true.
            return
          end if
        end do
      end do
    else
      do i = 1, size(b)
        do j = 1, size(a)
          if ( a(j) == b(i) ) then
            intersect = .true.
            return
          end if
        end do
      end do
    end if
    intersect = .false.

  end function Intersect

  ! -----------------------------------------------  Intersection  -----
  function Intersection ( A, B ) result ( C )
  ! Compute the intersection C of the sets A and B, each represented by
  ! arrays of integers.

    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error
    use Sort_M, only: Sort

    integer, intent(in) :: A(:), B(:)
    integer, pointer :: C(:) ! Intent(out) -- nullified and then allocated here

    integer :: I, J, K, Stat, TA(size(a)), TB(size(b)), TC(size(a)+size(b))

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
      MLSMSG_Allocate // 'C in Intersection' )
    c = tc(:k)

  end function Intersection

  ! ------------------------------------------------------  Union  -----
  function Union ( A, B ) result ( C )
  ! Compute the union C of the sets A and B, each represented by
  ! arrays of integers.

    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error
    use Sort_M, only: Sort

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

  end function Union

  ! --------------------------------------------------  UnionSize  -----
  integer function UnionSize ( A, B )
  ! Compute the size of the union of the sets A and B, each represented by
  ! arrays of integers.

    use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error
    use Sort_M, only: Sort

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

    unionSize = i

  end function UnionSize

! =====     Private Procedures     =====================================

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSSets

! $Log$
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
