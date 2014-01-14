! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module FIRST_SETS

  implicit NONE
  private

  integer, public, allocatable, save :: FIRST_PT(:) ! pointer to first set
  public :: FIND_FIRST

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine FIND_FIRST            ! find FIRST sets for nonterminal symbols

    use Error_m, only: Error
    use LISTS, only: ADDLTL, LCOMPR, LIST, NEW, REL
    use NULLABLE_M, only: NULLABLE
    use Tables, only: NTERMS, NUMPRD, NVOC, PRDIND => Prod_Ind, &
      & PRODCN => Productions

    logical :: ADD_CHANGE          ! "ADDLTL changed a list"
    logical :: CHANGE              ! "a first set changed"
    integer :: I, J, K             ! subscripts and loop inductors
    logical :: RESET(NTERMS+1:NVOC) ! "FIRST_PT(i) has been destroyed"

    allocate ( first_pt (nvoc) , stat=i )
    if ( i /= 0 ) call error &
      ( 'FIND_FIRST-E- Unable to allocate FIRST_PT, STAT =', 2 )

    ! Compute the first sets
    first_pt(nterms+1:nvoc) = 0    ! nonterminal first sets start empty
    do i = 1, nterms
      call new ( first_pt(i) )
      list(first_pt(i))%item = i        ! terminal first set is the terminal
    end do

    change = .true.
    do while ( change )
      change = .false.
      do i = 1, numprd             ! all the productions
        j = prdind(i)              ! beginning of i'th production
        do k = j+1, prdind(i+1)-1  ! all of i'th production's RHS
          if ( prodcn(k) /= prodcn(j) ) then ! j'th RHS symbol /= LHS
            call addltl ( first_pt(prodcn(k)), first_pt(prodcn(j)), add_change )
            if ( add_change ) change = .true.
          end if
          if ( .not. nullable(prodcn(k)) ) exit
        end do
      end do
    end do

    ! Now pack up the sets.  If two sets are identical they will be
    ! shared.  This is OK because the THEAD sets never change once
    ! generated.

    reset = first_pt(nterms+1:nvoc) == 0
    do i = nterms+1, nvoc-1
      if ( .not. reset(i) ) then
        do j = i+1, nvoc
          if ( .not. reset(j) ) then
            if ( lcompr(first_pt(i), first_pt(j)) ) then ! Sets equal?
              call rel (first_pt(j))    ! release one
              first_pt(j) = first_pt(i) ! use the other one
              reset(j) = .true.         ! indicate released one destroyed
            end if
          end if
        end do
      end if
    end do

    call print_first

  end subroutine FIND_FIRST

  subroutine PRINT_FIRST

    use Output_m, only: OUTPUT
    use Print_List, only: PNTLST
    use Processor_Dependent, only: NewPage
    use String_Table, only: Display_String, String_Length
    use Tables, only: First_Nonterminal, Last_Nonterminal, Vocab
    use TOGGLES_LR

    integer :: I, K                ! subscripts, loop inductors

    if (toggle(iachar('2')) /= 0) then
    ! Print the FIRST sets.
      call output ( newPage, dont_asciify=.true. )
      call output ( 'FIRST SETS:', advance='yes' )
      do k = First_Nonterminal, Last_Nonterminal
        if (first_pt(k) /= 0) then
          call display_string ( vocab(k), before='    ' )
          call output ( ': ' )
          i = string_length(vocab(k)) + 6
          call pntlst (first_pt(k), i, 10)
        end if
      end do
    end if

  end subroutine PRINT_FIRST

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module FIRST_SETS

! $Log$
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
