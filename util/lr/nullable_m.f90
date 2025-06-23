! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module NULLABLE_M

  implicit NONE
  private

  logical, public, allocatable, save :: NULLABLE(:)
  public :: FIND_NULLABLE

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine FIND_NULLABLE    ! find nullable nonterminals

    use Error_m, only: Error
    use Tables, only: First_Nonterminal, FRSPRD => First_Production, &
      & Last_Nonterminal, PRDIND => Prod_Ind, PRODCN => Productions
    use Toggles, only: GEN, Toggle
    use Trace, only: Trace_Begin, Trace_End

    logical :: CHANGE         ! "Something changed in NULLABLE"
    integer :: I, K           ! Subscripts, loop inductors

    if ( toggle(gen) ) call trace_begin ( 'Find_Nullable' )
    allocate ( nullable(last_nonterminal), stat=i )
    if ( i /= 0 ) call error &
      ( 'FIND_NULLABLE-E- Unable to allocate NULLABLE, STAT =', 2 )
    nullable = .false.        ! nothing nullable yet, esp. terminals
    do ! until .not. ichange
      change = .false.
      do k = First_Nonterminal, Last_Nonterminal   ! all the nonterminals
        if ( .not. nullable(k) ) then ! k is not yet known to be nullable
          i = frsprd(k)       ! First production with LHS = k
          if ( i > 0 ) then   ! Otherwise it an unused symbol defined by &V
            do while ( prodcn(prdind(i)) == k )  ! While LHS = k
              nullable(k) = all(nullable(prodcn(prdind(i)+1:prdind(i+1)-1)))
              if ( nullable(k) ) then
                change = .true.
                exit            ! go on to next nonterminal
              end if
              i = i + 1         ! go on to next RHS
            end do              ! While LHS = k
          end if
        end if
      end do                  ! all the nonterminals
      if ( .not. change ) exit
    end do
    call print_nullable
    if ( toggle(gen) ) call trace_end ( 'Find_Nullable' )

  end subroutine FIND_NULLABLE

  subroutine PRINT_NULLABLE
    use Output_m, only: Blanks, NewLine, OUTPUT
    use Processor_Dependent, only: NewPage
    use String_Table, only: Display_String, String_Length
    use Tables, only: First_Nonterminal, Last_Nonterminal, VOCAB
    use TOGGLES_LR
    integer :: I, K, L           ! Subscripts, loop inductors
    if (toggle(iachar('2')) /= 0) then
    ! Print the list of nullable nonterminals.
      call output ( newPage, dont_asciify=.true. )
      call output ( 'POTENTIALLY NULL NON-TERMINALS: ', advance='yes' )
      call blanks ( 4 )
      i = 4
      do k = First_Nonterminal, Last_Nonterminal   ! all the nonterminals
        if ( nullable(k) ) then
          l = string_length(vocab(k))
          if ( l+i > 120 ) then
            call newLine
            call blanks ( 4 )
          end if
          call display_string ( vocab(k), before=' ' )
          i = i + 1 + l
        end if
      end do
      call newLine
    end if

  end subroutine PRINT_NULLABLE

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module NULLABLE_M

! $Log$
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
