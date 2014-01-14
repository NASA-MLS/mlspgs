! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Print_The_Grammar_m

  ! Flatten the abstract syntax tree for a grammar to the data structures
  ! used to analyze the grammar.

  implicit NONE
  private

  public :: Print_The_Grammar

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Procedures     =====================================

  subroutine Print_The_Grammar ( Prod_Ind, Productions, Actions, P )
  ! Pretty-print the Production and Actions tables
    use Processor_Dependent, only: NewPage
    use String_Table, only: Display_String, String_Length
    use Output_m, only: Blanks, Newline, Output
    integer, intent(in) :: Prod_Ind(:) ! Production I is in Productions
                             ! Productions(Prod_Ind(i):Prod_Ind(i+1)-1),
                             ! including the LHS in every one.
    integer, intent(in) :: Productions(:)
    integer, intent(in) :: Actions(:)  ! Actions(i) is the action for the
                             ! I'th production
    integer, intent(in) :: P(:)        ! String index of symbol with index
                             ! given in the Productions array

    integer :: I, J          ! Loop inductors
    integer :: LHS, Prev_LHS ! Left-hand-side string index
    integer :: W             ! Width of left-hand-side printing area

    call output ( newPage // '     T H E   P R O D U C T I O N S', advance='yes', &
      & dont_asciify = .true. )
    prev_lhs = 0
    do i = 1, size(prod_ind,1)-1
      lhs = p(productions(prod_ind(i)))
      if ( lhs /= prev_lhs ) call newLine
      call output ( i, 4 )
      if ( lhs /= prev_lhs ) then
        call display_string ( lhs, before=' ' )
        w = 1 + string_length(lhs )
        prev_lhs = lhs
      else
        call blanks ( w )
      end if
      call output ( ' ->' )
      do j = prod_ind(i)+1, prod_ind(i+1)-1
        if ( productions(j) == 0 ) then
          call output ( ' <my goal>' )
        else
          call display_string ( p(productions(j)), before=' ' )
        end if
      end do
      if ( actions(i) /= 0 ) then
        call output ( ' =>' )
        call display_string ( abs(actions(i)), before=' ' )
        if ( actions(i) < 0 ) call output ( ' ?' )
      end if
      call newLine
    end do

  end subroutine Print_The_Grammar

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Print_The_Grammar_m

! $Log$
