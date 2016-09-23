! Copyright 2015, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module QTM_dumps
!=============================================================================

  ! Dump routines to go with QTM_m.  They're here instead of there so that
  ! the output module and all its support aren't dragged into QTM_m.

  implicit NONE
  private

  public :: Dump, Dump_QID, Dump_Stack

  interface Dump
    module procedure Dump_Stack
  end interface

  character(len=*), parameter, public :: StackDumpHead = &
    & ' QL  pNode X         Y      xNode X         Y      yNode X         Y       QID'

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Dump_QID ( QID, Before, After, Advance )
    use Output_m, only: Output
    use QTM_m, only: High_Bit_Index, QK
    integer(qk), intent(in) :: QID
    character(*), intent(in), optional :: Before, After, Advance
    integer :: H
    if ( present(before) ) call output ( before )
    if ( qid < 8 ) then
      call output ( qid, format='(z0)' )
    else
      h = high_bit_index ( QID ) - 4
      call output ( int(shiftr(qid,h)) - 7 ) ! Octant number
      do while ( h /= 0 )
        h = h - 2
        call output ( iand(int(shiftr(qid,h)),3) )
      end do
      if ( present(after) ) call output ( after )
    end if
    call output ( "", advance )
  end subroutine Dump_QID

  subroutine Dump_Stack ( Stack, Oct, Top, Advance, Head )
    use Output_m, only: Output
    use QTM_m, only: Stack_t
    type(stack_t), intent(in) :: Stack
    integer, intent(in), optional :: Oct
    logical, intent(in), optional :: Top  ! Dump only the top frame
    character(*), intent(in), optional :: Advance
    character(*), intent(in), optional :: Head ! Print the heading followed
                                               ! by trim(head)
    integer :: I, I1, L, O, Oct1, Octn, Pn, Xn, Yn
    character(127) :: Line
    character(3) :: MyAdv
    logical :: MyTop
    oct1 = 8
    octn = 15
    if ( present(oct) ) then
      oct1 = oct
      octn = oct
    end if
    myTop = .false. ! Dump the whole stack
    if ( present(top) ) myTop = top
    do o = oct1, octn
      i1 = 1
      if ( myTop ) i1 = stack%top(o)
      if ( present(head) ) then
        call output ( stackDumpHead )
        call output ( trim(head), advance=advance )
      end if
      do i = i1, stack%top(o)
        xn = stack%xNode(i,o)
        yn = stack%yNode(i,o)
        pn = 6 - ( xn + yn )
        write ( line, '(i3,4(i3,2f10.5),i3)' ) i, &
          & pn, stack%z(pn,i,o), &
          & xn, stack%z(xn,i,o), &
          & yn, stack%z(yn,i,o)
        l = len_trim(line) + 4
        if ( myTop ) then
          write ( line(l:), '(99i0)' ) stack%qid(1,o)-7, &
            & stack%qid(2:stack%top(o),o)
        else
          write ( line(l:), '(i0)' ) stack%qid(i,o) - merge(7,0,i==1)
        end if
        myAdv = 'yes'
        if ( i == stack%top(o) ) myAdv = advance
        call output ( trim(line), myAdv )
      end do
    end do
  end subroutine Dump_Stack

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module QTM_dumps

! $Log$
! Revision 2.1  2016/09/23 01:54:18  vsnyder
! Initial commit; stuff moved here from QTM_m
!
