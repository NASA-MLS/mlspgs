! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Xref

  ! Flatten the abstract syntax tree for a grammar to the data structures
  ! used to analyze the grammar.

  implicit NONE
  private

  public :: Cross_Reference

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Procedures     =====================================

  subroutine Cross_Reference ( Productions, Prod_Ind, Strings )

  ! Print a cross-reference of vocabulary symbols.  Positive integer is
  ! reference in RHS, negative is reference as LHS.

    use Output_m, only: Blanks, NewLine, Output
    use Tables, only: First_Nonterminal, First_Terminal, Last_Nonterminal
    use Processor_Dependent, only: NewPage
    use String_Table, only: Display_String, String_Length

    integer, intent(in) :: Productions(:)
    integer, intent(in) :: Prod_Ind(:) ! Production I is in
                                       ! productions(prod_ind(i):prod_ind(i+1)-1)
    integer, intent(in) :: Strings(:)  ! Mapping from sorted to original string indices


    type :: Ref
      integer :: Prod = 0
      type(ref), pointer :: Prev => NULL()
    end type Ref

    integer :: I, J, L, LStart
    type(ref), pointer :: NewRef
    type(ref) :: Refs(first_terminal:last_nonterminal)

    ! Create lists of symbol references
    do i = 1, size(prod_ind)-1
      allocate ( newRef )
      newRef%prod = -i
      j = prod_ind(i)
      newRef%prev => refs(productions(j))%prev
      refs(productions(j))%prev => newRef
      do j = prod_ind(i)+1, prod_ind(i+1)-1
        allocate ( newRef )
        newRef%prod = i
        newRef%prev => refs(productions(j))%prev
        refs(productions(j))%prev => newRef
      end do
    end do

    ! Print the cross reference
    call output ( newPage, dont_asciify=.true. )
    call blanks ( 30 )
    call output ( 'A   V O C A B U L A R Y   C R O S S   R E F E R E N C E', &
      & advance='yes' )

    do i = first_terminal, last_nonterminal
      if ( i == first_terminal ) then
        call output ( ' T E R M I N A L S', advance='yes' )
        call newLine
      else if ( i == first_nonterminal ) then
        call newLine
        call output ( ' N O N T E R M I N A L S', advance='yes' )
        call newLine
      end if
      call display_string ( strings(i), before=' ' )
      l = string_length ( strings(i) )
      lstart = l+5 - mod(l+5,5)
      call blanks ( lstart-l )
      l = lstart
      if ( associated(refs(i)%prev) ) then
        call print_the_refs ( refs(i)%prev )
      else
        call output ( '    Is not referenced; this is probably an error', advance='yes' )
      end if
      call newLine
    end do

  contains

    recursive subroutine Print_The_Refs ( Item )
    ! Print the references linked to Ref.  The list is in reverse order,
    ! so go to the tail of it first, and then print and deallocate as
    ! the recursion unwinds.

      type(ref), pointer :: Item
      character(5) :: Text

      if ( associated(item%prev) ) call print_the_refs ( item%prev )
      if ( item%prod /= 0 ) then
        write ( text, '(i5)' ) item%prod
        l = l + 5
        if ( l > 115 ) then
          call newLine
          call blanks ( lstart+1 )
          l = lstart + 1
        end if
        call output ( text )
        deallocate ( item )
      end if

    end subroutine Print_The_Refs

  end subroutine Cross_Reference

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Xref

! $Log$
! Revision 1.1  2014/01/14 00:15:14  vsnyder
! Initial commit of new module for new LR
!
