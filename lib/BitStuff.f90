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
module BitStuff
!=============================================================================

! This module contains routines for interrogating, manipulating, using bits.

  use MLSSets, only: FindAll

  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

! In the following "bits" refers in most cases to the bits
! held within an integer. Thus an integer can be decomposed
! into an array of bits. How large an array and what the starting index to use
! are very important. We will assume that the starting index is 0 and 
! that the largest index is MAXBITNUMBER.

!     (subroutines and functions)
! BitsToBooleans     Converts bits to an array of Boolean
! BooleansToBits     Converts an array of Boolean to bits
! CountBits          Returns the number of non-zero bits
! IsBitPatternSet    Are any of the bits set in bitpattern also set in the arg
! IsBitSet           Is the specified bit number set
! NegativeIfBitPatternSet  Return the negative abs. of the first arg
!                      if any of the bits set in bitpattern also set in i
! NegativeIfBitSet   Return the negative abs. of the first arg
!                      if the bitnumber is set in i
! SetBits            Return the integer represented by the bits set among set
! WhichBitsAreSet    Which bit numbers are set?
! === (end of toc) ===
! === (start of api) ===
! BitsToBooleans ( bits i, log Booleans(:) )
! BooleansToBits ( log Booleans(:), bits i )
! int CountBits ( value )
!                 value can be a scalar or an array of ints or chars
! log IsBitPatternSet ( bits i, bits bitPattern )
! log IsBitSet ( bits i, int bitNumber )
! value NegativeIfBitPatternSet ( value, bits i, bits bitPattern )
! value NegativeIfBitSet ( value, bits i, int bitNumber )
!                  value can be an int, real or double, scalar or array
! bits SetBits ( int set(:) )
! WhichBitsAreSet ( bits i, int set(:), [int howMany], [int unset(:)] )
! === (end of api) ===

  public :: BitsToBooleans, BooleansToBits
  public :: CountBits, CountBits_0, CountBits_1, CountBits_2
  public :: CountCharBits_0, CountCharBits_1, CountCharBits_2
  public :: isBitSet, IsBitPatternSet, SetBits, WhichBitsAreSet
  public :: NegativeIfBitSet, NegativeIfBitPatternSet
  public :: Reverse

  interface CountBits
    module procedure countBits_0, countBits_1, countBits_2
    module procedure countCharBits_0, countCharBits_1, countCharBits_2
  end interface

  interface NegativeIfBitSet
    module procedure NegativeIfBitSet_int, NegativeIfBitSet_sng, &
      & NegativeIfBitSet_dbl
  end interface

  interface NegativeIfBitPatternSet
    module procedure NegativeIfBitPatternSet_int, NegativeIfBitPatternSet_sng, &
      & NegativeIfBitPatternSet_dbl
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! bitNumber starts with 0 at the smallest
  ! and goes up from there
  integer, private :: ibitindx ! A dummy index used in array constructor
  
  ! Note: the following is assumed to be the max bit number by all
  ! procedures in this module except for Reverse
  ! which relies on the intrinsic bit_size which should return
  ! an integer = (MAXBITNUMBER+1)
  ! If you plan on using the software on a platform for which
  ! bit_size returns a different value, change the following
  integer, parameter, public :: MAXBITNUMBER = 30 ! The 31st is the sign bit
  integer, parameter, dimension(0:MAXBITNUMBER) :: BITINDX = &
    & (/ (ibitindx, ibitindx = 0, MAXBITNUMBER) /)

contains

  ! ------------------------------------------------  BitsToBooleans  -----
  subroutine BitsToBooleans ( i, Booleans )
    ! Convert an integer to an array of booleans, treating
    ! the integer as a bit pattern
    ! E.g., given 41 (101001 in binary) returns (/ T, F, F, T, F, T /)
    ! Note that the first element of the array Booleans
    ! corresponds to bit number 0
    ! That's why we begin its index number with 0 instead of 1
    integer, intent(in)                 :: i
    logical, dimension(0:), intent(out) :: Booleans
    ! Local variables
    integer :: bitNumber
    ! Executable
    Booleans = .false.
    do bitNumber = 0, min(size(Booleans) - 1, MAXBITNUMBER)
      Booleans(bitNumber) = isBitSet( i, bitNumber )
    enddo
  end subroutine BitsToBooleans

  ! ------------------------------------------------  BooleansToBits  -----
  subroutine BooleansToBits ( Booleans, i )
    ! Convert an an array of booleans to an integer, treating
    ! the integer as a bit pattern
    ! E.g., given (/ T, F, F, T, F, T /) returns 41 (101001 in binary) 
    ! Note that the first element of the array Booleans
    ! corresponds to bit number 0
    ! That's why we begin its index number with 0 instead of 1
    logical, dimension(0:), intent(in)  :: Booleans
    integer, intent(out)                :: i
    ! Local variables
    integer :: bitNumber
    ! Executable
    i = 0
    do bitNumber = 0, min( ubound(Booleans, 1), bit_size(i) - 1 )
      if ( Booleans(bitNumber) ) i = ibset( i, bitNumber )
    enddo
  end subroutine BooleansToBits

  ! ------------------------------------------------  CountBits_0  -----
  integer function CountBits_0 ( Word )
    integer, intent(in) :: Word
    ! Count the number of nonzero bits in Word.
    integer :: Counts(0:255) = & ! Number of nonzero bits in each byte
      !    0  1  2  3  4  5  6  7  8  9  A  B  C  D  E  F
      & (/ 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,  & ! 0 
      &    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,  & ! 1 
      &    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,  & ! 2 
      &    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,  & ! 3 
      &    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,  & ! 4 
      &    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,  & ! 5 
      &    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,  & ! 6 
      &    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,  & ! 7 
      &    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,  & ! 8 
      &    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,  & ! 9 
      &    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,  & ! A 
      &    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,  & ! B 
      &    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,  & ! C 
      &    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,  & ! D 
      &    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,  & ! E 
      &    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8 /)  ! F 
    integer :: I

    countBits_0 = counts(iand(word,255))
    do i = 8, bit_size(word)-8, 8
      countBits_0 = countBits_0 + counts(iand(ishft(word,-i),255))
    end do
  end function CountBits_0

  ! ------------------------------------------------  CountBits_1  -----
  integer function CountBits_1 ( Array )
    integer, intent(in) :: Array(:)
    ! Count the number of nonzero bits in all of the elements of Word.
    integer :: I
    countBits_1 = countBits(array(1))
    do i = 2, size(array)
      countBits_1 = countBits_1 + countBits(array(i))
    end do
  end function CountBits_1

  ! ------------------------------------------------  CountBits_2  -----
  integer function CountBits_2 ( Array )
    integer, intent(in) :: Array(:,:)
    ! Count the number of nonzero bits in all of the elements of Word.
    integer :: I
    countBits_2 = countBits(array(:,1))
    do i = 2, size(array,2)
      countBits_2 = countBits_2 + countBits(array(:,i))
    end do
  end function CountBits_2

  ! --------------------------------------------  CountCharBits_0  -----
  integer function CountCharBits_0 ( Char, What )
    ! Count the number of nonzero bits.
    character(len=1), intent(in) :: Char
    integer, intent(in), optional :: What ! mask of which bits to count
    integer :: myWhat
    myWhat = not(0)
    if ( present(what) ) myWhat = what
    countCharBits_0 = countBits(iand(ichar(char),myWhat))
  end function CountCharBits_0

  ! --------------------------------------------  CountCharBits_1  -----
  integer function CountCharBits_1 ( Array, What )
    character(len=1), intent(in) :: Array(:)
    integer, intent(in), optional :: What ! mask of which bits to count
    ! Count the number of nonzero bits in all of the elements of Word.
    integer :: I
    countCharBits_1 = countBits(array(1),what)
    do i = 2, size(array)
      countCharBits_1 = countCharBits_1 + countBits(array(i),what)
    end do
  end function CountCharBits_1

  ! --------------------------------------------  CountCharBits_2  -----
  integer function CountCharBits_2 ( Array, What )
    character, intent(in) :: Array(:,:)
    integer, intent(in), optional :: What ! mask of which bits to count
    ! Count the number of nonzero bits in all of the elements of Word.
    integer :: I
    countCharBits_2 = countBits(array(:,1),what)
    do i = 2, size(array,2)
      countCharBits_2 = countCharBits_2 + countBits(array(:,i),what)
    end do
  end function CountCharBits_2

  ! --------------------------------------------  IsBitPatternSet  -----
  elemental function IsBitPatternSet ( i, bitPattern ) result(SooDesu)
    ! Are any of the bitNumbers in bitPattern set in i?
    ! E.g., if bitPattern=1 bitNumber is 0 and we
    ! simply test whether i is odd
    integer, intent(in) :: i
    integer, intent(in) :: bitPattern
    logical             :: SooDesu
    ! Executable
    SooDesu = ( iand(i, bitPattern) /= 0 )
  end function IsBitPatternSet

  ! --------------------------------------------  IsBitSet  -----
  elemental function IsBitSet ( i, bitNumber ) result(SooDesu)
    ! Is the bitNumber bit of i set?
    ! if omitted, bitNumber is assumed 0 and we
    ! simply test whether i is odd
    integer, intent(in) :: i
    integer, intent(in), optional :: bitNumber
    logical             :: SooDesu
    integer :: myBit
    ! Executable
    myBit = 0
    if ( present(bitNumber) ) myBit = bitNumber
    SooDesu = btest(i, myBit)
  end function IsBitSet

  ! --------------------------------------------  NegativeIfBitPatternSet  -----
  ! This family of functions returns the negative abs. val.
  ! of value if any of the bitNumbers in bitPattern set in i
  ! This may be useful if precisions for certain radiances are to be
  ! set negative (which means don't use them) when a certain bright object
  ! is in the field of view
  elemental function NegativeIfBitPatternSet_int ( value, i, bitPattern ) &
    & result(res)
    integer, intent(in) :: value
    integer, intent(in) :: i
    integer, intent(in) :: bitPattern
    integer             :: res
    ! Executable
    res = value
    if ( IsBitPatternSet(i, bitPattern) ) res = - abs(value)
  end function NegativeIfBitPatternSet_int

  elemental function NegativeIfBitPatternSet_sng ( value, i, bitPattern ) &
    & result(res)
    real, intent(in)    :: value
    integer, intent(in) :: i
    integer, intent(in) :: bitPattern
    real                :: res
    ! Executable
    res = value
    if ( IsBitPatternSet(i, bitPattern) ) res = - abs(value)
  end function NegativeIfBitPatternSet_sng

  elemental function NegativeIfBitPatternSet_dbl ( value, i, bitPattern ) &
    & result(res)
    double precision, intent(in) :: value
    integer, intent(in) :: i
    integer, intent(in) :: bitPattern
    double precision    :: res
    ! Executable
    res = value
    if ( IsBitPatternSet(i, bitPattern) ) res = - abs(value)
  end function NegativeIfBitPatternSet_dbl

  ! --------------------------------------------  NegativeIfBitSet  -----
  ! This family of functions returns the negative abs. val.
  ! of value if the bitNumber of i is set
  ! This may be useful if precisions for certain radiances are to be
  ! set negative (which means don't use them) when a certain bright object
  ! is in the field of view
  elemental function NegativeIfBitSet_int ( value, i, bitNumber ) result(res)
    integer, intent(in) :: value
    integer, intent(in) :: i
    integer, intent(in) :: bitNumber
    integer             :: res
    ! Executable
    res = value
    if ( isBitSet(i, bitNumber) ) res = - abs(value)
  end function NegativeIfBitSet_int

  elemental function NegativeIfBitSet_sng ( value, i, bitNumber ) result(res)
    real, intent(in)    :: value
    integer, intent(in) :: i
    integer, intent(in) :: bitNumber
    real                :: res
    ! Executable
    res = value
    if ( isBitSet(i, bitNumber) ) res = - abs(value)
  end function NegativeIfBitSet_sng

  elemental function NegativeIfBitSet_dbl ( value, i, bitNumber ) result(res)
    double precision, intent(in) :: value
    integer, intent(in) :: i
    integer, intent(in) :: bitNumber
    double precision    :: res
    ! Executable
    res = value
    if ( isBitSet(i, bitNumber) ) res = - abs(value)
  end function NegativeIfBitSet_dbl

  ! --------------------------------------------  SetBits  -----
  function Reverse ( i ) result(anti)
    ! Return the integer reversing the bit order of i
    ! E.g., if set(i) = (/ 0, 2, 4 /) set(result) will be (/N, N-2, N-4/)
    ! where N = MAXBITNUMBER(i)
    
    ! Because the 31st bit is the sign bit, the zeroth bit would be discarded
    ! Therefore we multiply by 2 before reversing
    
    ! Restrictions:
    ! i >= 0
    ! 2*i < huge(i)
    
    ! Always returns 0 if restrictions violated
    integer, intent(in) :: i
    integer :: anti
    ! Local variables
    integer :: howmany
    integer :: N
    integer, dimension(bit_size(1)) :: set
    ! Executable
    anti = 0
    ! Error checks
    if ( i < 0 .or. i > huge(i)/2 ) return
    N = bit_size(i) - 1 ! Because bit numbers start with 0, 31st is sign
    call WhichBitsAreSet ( 2*i, set, howmany )
    if ( howmany == 0 ) return
    set = N - set
    anti = SetBits(set)
  end function Reverse

  ! --------------------------------------------  SetBits  -----
  function SetBits ( set ) result(i)
    ! Return the integer represented by setting the bitNumbers
    ! in the array set
    ! E.g., if set = (/ 0, 2, 4 /) result will be 21 (= b'10101')
    integer, intent(in), dimension(:) :: set
    integer :: i
    integer :: indx
    integer :: iset
    ! Executable
    i = 0
    if ( size(set) < 1 ) return
    do iset = 1, size(set)
      indx = set(iset)
      if ( indx < 0 .or. indx > MAXBITNUMBER ) cycle
      i = ibset(i, indx)
    enddo
  end function SetBits

  ! --------------------------------------------  WhichBitsAreSet  -----
  subroutine WhichBitsAreSet ( i, set, howmany, unset )
    ! Which bitNumbers of i are set? Returned in array "set"
    ! E.g., given 41 (101001 in binary) returns (/ 0,3,5 /)
    ! If present, unset is filled with bit numbers not set (/ 1,2,4,6,7,.. /)
    ! Note: because 0 is an actual bit number, we prefill
    ! arrays with -1. So when returned array has count(set > -1)
    ! elements you will know how to proceed. (No, how?)
    
    ! NAG has caused us to employ unsightly array constructors
    ! because it failed to properly elementalize the logical function
    ! IsBitSet
    integer, intent(in)                          :: i
    integer, intent(out), dimension(:) :: set
    integer, intent(out), optional               :: howmany
    integer, intent(out), dimension(:), optional :: unset
    integer :: N, indx
    ! Executable
    if ( present(unset) ) then
      N = max( size(set), size(unset) )
      N = min( N, MAXBITNUMBER+1 )
      set = -1
      unset = -1
      call FindAll( &
        & (/ (IsBitSet(i, indx), indx=0, N) /), &
        & set, howmany, unset )
      ! The first bit number is 0, not 1
      set = set - 1
      unset = unset - 1
    else
      N = min( size(set), MAXBITNUMBER+1 )
      set = -1
      call FindAll( &
        & (/ (IsBitSet(i, indx), indx=0, N) /), &
        & set, howmany )
      set = set - 1
    endif
  end subroutine WhichBitsAreSet

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module BitStuff

! $Log$
! Revision 2.11  2007/01/31 00:06:23  pwagner
! Made MAXBITNUMBER public
!
! Revision 2.10  2006/02/06 22:40:11  pwagner
! Added Reverse function
!
! Revision 2.9  2006/02/01 23:41:05  pwagner
! Can convert bits <-> booleans
!
! Revision 2.8  2005/11/15 00:17:20  pwagner
! Changes to workaround NAG failure to elemntalize properly
!
! Revision 2.7  2005/11/11 21:42:17  pwagner
! Added bit setting, retrieval, and conditional sign-changing procedures
!
! Revision 2.6  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.5  2003/05/10 00:56:31  livesey
! Bug fix in CountCharBits_1/2
!
! Revision 2.4  2002/10/08 00:09:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.3  2002/02/13 20:52:26  vsnyder
! Added a 'what' mask to count_bits_char...
!
! Revision 2.2  2002/02/05 02:39:59  vsnyder
! Change mask from 1-bit per to 8-bits per (using character)
!
! Revision 2.1  2001/10/03 17:36:15  vsnyder
! Initial commit
!
