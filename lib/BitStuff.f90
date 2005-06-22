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

! This module contains routines for manipulating bits.

  implicit NONE
  private
  public :: CountBits, CountBits_0, CountBits_1, CountBits_2
  public :: CountCharBits_0, CountCharBits_1, CountCharBits_2

  interface CountBits
    module procedure countBits_0, countBits_1, countBits_2
    module procedure countCharBits_0, countCharBits_1, countCharBits_2
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
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
