! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module BitStuff
!=============================================================================

! This module contains routines for manipulating bits.

  implicit NONE
  private
  public :: CountBits, CountBits_0, CountBits_1, CountBits_2
  public :: CountCharBits_0, CountCharBits_1, CountCharBits_2

  interface CountBits
    module procedure countBits_0, countBits_1, countBits_2
    module procedure countCharBits_0, countCharBits_1, countCharBits_2
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

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
    countCharBits_1 = countBits(array(1))
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
    countCharBits_2 = countBits(array(:,1))
    do i = 2, size(array,2)
      countCharBits_2 = countCharBits_2 + countBits(array(:,i),what)
    end do
  end function CountCharBits_2

end module BitStuff

! $Log$
! Revision 2.3  2002/02/13 20:52:26  vsnyder
! Added a 'what' mask to count_bits_char...
!
! Revision 2.2  2002/02/05 02:39:59  vsnyder
! Change mask from 1-bit per to 8-bits per (using character)
!
! Revision 2.1  2001/10/03 17:36:15  vsnyder
! Initial commit
!
