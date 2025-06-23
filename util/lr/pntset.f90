! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Print_Set

  implicit NONE
  private
  public :: PNTSET

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine PNTSET ( LNADQT )

    use Basis_m, only: BASIS, IFINAL, INDBAS, Items, Item_t
    use Complete, only: Closures_Index, Closures
    use Delete, only: DELCS
    use Error_m, only: Error
    use Output_m, only: Blanks, NewLine, OUTPUT
    use Processor_Dependent, only: NewPage
    use LISTS, only: LINT, LIST, NEW, REL
    use Print_List, only: PNTLST
    use String_Table, only: Display_String, String_Length
    use Tables, only: Actions, NUMPRD, PRDIND => Prod_Ind, &
      & PRODCN => Productions, VOCAB
    use TOGGLES_LR
    use Transitions_And_Reductions, only: Red, Tran

    integer, intent(out) :: LNADQT     ! The last inadequate state.  Zero if
                             ! the grammar is adequate, i.e., it is LR(1).

    ! Print the configuration sets.

    ! *****     External References     ************************************

    ! DELCS   deletes a reference to a context set.
    ! ERROR   prints error messages.
    ! LINT    detects intersections between pairs of lists.
    ! PNTLST  prints a list of vocabulary items.

    ! *****     Local Variables     ****************************************

    ! ADEQUT  indicates whether the current state is adequate, that is, free
    !         of conflicts.
    ! I       points to the current basis.
    ! IBASE   is the base of the set of productions being printed.
    ! IPR     points to the set of productions for a configuration set.
    ! IPTR    points to a context set list element.
    ! ITEMP   points to a temporary list node.
    ! J       is a loop inductor and subscript for RED.
    ! JDOT    is the position being printed in the production.
    ! JEND    is the upper limit for J when reductions are being printed, or
    !         the upper limit for JSTART when transitions are being printed.
    ! JSTART  is the starting value for J when reductions are being printed,
    !         or the starting value, loop inductor, and subscript for TRAN
    !         when transitions are being printed.
    ! K       is a loop induction variable and subscript.
    ! L       is a loop inductor, subscript for RED when conflicts are
    !         being detected, and position in output line.
    ! LINK    is a pointer to a list node when cross reference lists are
    !         being constructed or printed.
    ! LSTHED  is an array of pointers to lists of states in which a
    !         a production is analyzed.
    ! REDUCE  is .TRUE. iff the last configuration printed was a reducing
    !         configuration, in which case the context is not printed (since
    !         it is the same as the lookahead).

    logical ADEQUT
    integer I, IBASE, IPR, IPTR, ITEMP, JDOT, JEND, JSTART
    integer K, L, LINK
    integer LSTHED(NUMPRD)
    logical REDUCE

  !     *****     Procedures     *****************************************

    if (toggle(iachar('M')) /= 0) then
      call clean_up ! redundant context set references
      return
    end if
    lnadqt = 0
    call new (itemp)
    call output ( newPage, dont_asciify=.true. )
    call blanks ( 37 )
    call output ( 'T H E   P A R S I N G   A U T O M A T O N', advance='yes' )
    call newLine
    call output ( ' State number', advance='yes' )
    call output ( '     .  Production numbers', advance='yes' )
    call output ( '          .  State Contents', advance='yes' )
    do i = 1, indbas-1
      ipr = basis(i)%item        ! first item for the state
      call print_the_basis       ! for the state
      if ( toggle(iachar('C')) /= 0 ) call print_the_closure
      call print_the_transitions ! from the state
      call print_the_reductions  ! from the state
    end do
    call rel (itemp)
    if (lnadqt /= 0) then
      call error ( 'This grammar is not LR(1)', 1 )
      call output ( lnadqt, before=' Last inadequate state:', advance='yes' )
    end if
    if (toggle(iachar('X')) == 0) &
      call print_cross_reference ! of productions and states

  contains

    subroutine clean_up ! redundant context set references
      integer :: I, J
      do i = 1, indbas-1
        ipr = basis(i)%item
        if ( items(ipr)%prod == 1 .and. items(ipr)%dot > 3 ) ifinal = i
        do j = ipr, basis(i+1)%item-1
          call delcs (items(j)%set)
        end do
      end do
    end subroutine clean_up ! redundant context set references

    subroutine print_the_basis ! for the state
      integer :: J
      integer :: L ! How many blanks before printing the production.
      call newLine
      call output ( i, 6 )
      if ( items(ipr)%prod == 1 .and. items(ipr)%dot > 3 ) ifinal = i
      l = 0   ! State number has been printed; don't need initial space
      do j = ipr, basis(i+1)%item - 1
        call print_the_production ( items(j), l )! with a dot in the right side
        l = 6 ! No state number after additional productions; need initial space
        if (.not.reduce) then
          iptr = list(items(j)%set)%next
          if (iptr /= 0 .and. toggle(iachar('2')) /= 0) then
            call blanks ( 12 )
            call output ( 'Context:' )
            call pntlst (iptr, 23, 25)
          end if
        end if
        call delcs ( items(j)%set )
      end do
    end subroutine print_the_basis ! for the state

    subroutine print_the_closure ! for the state
      integer :: J      ! Closure to print
      integer :: J1, J2 ! Range of Closures to print
      integer :: NBasis ! Number of items in the basis
      nBasis = basis(i+1)%item - basis(i)%item
      j1 = closures_index(i-1) + nBasis + 1
      j2 = closures_index(i)
      if ( j1 > j2 ) return
      call blanks ( 12 )
      call output ( 'Closure:', advance='yes' )
      do j = j1, j2
        call print_the_production ( closures(j), 6 )
      end do
    end subroutine print_the_closure

    subroutine print_the_production ( Item, Initial )! with a dot in the right side

      type(item_t), intent(in) :: Item  ! Configuration set item
      integer, intent(in) :: Initial    ! Space before the production

      integer :: K
      integer :: Prod_Number
      integer :: W ! String length of vocab item

      ! Print the left side.

      prod_number = item%prod
      ibase = prdind(prod_number)
      call blanks ( initial )
      call output ( prod_number, 5 )
      call display_string ( vocab(prodcn(ibase)), before=' ' )
      call output ( ' ->' )
      l = string_length(vocab(prodcn(ibase))) + 16

      ! Print the right side of the production with a dot before the
      ! IDOT'th item.  IDOT is in item%dot.

      jdot = 1
      do k = ibase+1, prdind(prod_number+1)-1
        w = string_length(vocab(prodcn(k)))
        if ( w + l > 120 ) then
          call newLine
          call blanks ( 13 )
          l = 13
        end if
        if ( jdot == item%dot ) then
          call output ( ' .' )
          l = l + 2
        end if
        call display_string ( vocab(prodcn(k)), before=' ' )
        l = l + w + 1
        jdot = jdot + 1
      end do
      reduce = jdot == item%dot
      if (reduce) then
        call output ( ' .' )
        if ( actions(prod_number) /= 0 ) then
          call display_string ( abs(actions(prod_number)), before=' => ' )
          if ( actions(prod_number) < 0 ) call output ( ' ?' )
        end if
      end if
      call newLine
    end subroutine print_the_production ! with a dot in the right side

    subroutine print_the_transitions ! from the state
      jstart = basis(i)%tran
      jend = basis(i+1)%tran - 1
      if (jstart <= jend) then
        call blanks ( 12 )
        call output ( 'Transitions: ' )
        l = 12 + len('Transitions: ')
        do while (jstart <= jend)
          if (l > 115) then
            call newLine
            call blanks ( 16 )
          end if
          call output ( tran(jstart), 5 )
          l = l + 5
          jstart = jstart + 1
        end do
        call newLine
      end if
    end subroutine print_the_transitions ! from the state

    subroutine print_the_reductions ! from the state
      integer :: J
      adequt = .true.
      jstart = basis(i)%red
      jend = basis(i+1)%red - 1
      if (jstart <= jend) then
        call blanks ( 12 )
        call output ( 'Reductions:' )
        do j = jstart, jend
          if ( j /= jstart ) call blanks ( 12+len('Reductions:'))
          call output ( red(j)%prod, 7 )
          iptr = list(red(j)%set)%next
          call pntlst (iptr, 30, 30)
          call check_for_conflicts ( j ) ! in state j
        end do
        if (.not. adequt) then
          if (lnadqt /= 0) then
            call output ( lnadqt, before=' Last inadequate state:', advance='yes' )
          end if
          lnadqt = (i+4)/5
        end if
      end if
    end subroutine print_the_reductions ! from the state

    subroutine check_for_conflicts ( j ) ! in state j
      integer, intent(in) :: J
      integer :: L
      do l = basis(i)%tran, basis(i+1)%tran-1
        k = basis(tran(l))%item
        list(itemp)%item = prodcn(prdind(items(k)%prod)+items(k)%dot-1)
        if ( lint(itemp, list(red(j)%set)%next) ) then
          adequt = .false.
          call output ( (tran(l)+4)/5, &
            & before=' *** Intersection with transition to ' )
          call display_string ( vocab(list(itemp)%item), before=' on ', advance='yes' )
        end if
      end do
      if (j /= jstart) then
        do l = jstart, j-1
          if ( lint(list(red(l)%set)%next, list(red(j)%set)%next) ) then
            adequt = .false.
            call output ( red(l)%prod, before=' *** Intersection with reduction of ', &
              & advance='yes' )
          end if
        end do
      end if
    end subroutine check_for_conflicts ! in this state

    subroutine print_cross_reference ! of productions and states
      integer :: I, J
      lsthed(1:numprd) = 0
      do i = 1, indbas-1
        do j = basis(i)%item, basis(i+1)%item - 1
          call new (link)
          list(link)%item = i
          list(link)%next = lsthed(items(j)%prod)
          lsthed(items(j)%prod) = link
        end do
      end do

      ! Print the cross reference lists.

      call output ( newPage, dont_asciify=.true. )
      call blanks ( 34 )
      call output ( 'A U T O M A T O N   C R O S S   R E F E R E N C E', &
        & advance='yes' )
      call newLine
      call output ( 'Production', advance='yes' )
      call output ( '    .    Analyzed in states', advance='yes' )
      do i = 1, numprd
        link = lsthed(i)
        if (link /= 0) then
          call output ( i, 5 )
          l = 6
          do while (link /= 0)
            if (l > 115) call newLine
            call output ( list(link)%item, 5 )
            l = l + 5
            link = list(link)%next
          end do
          call newLine
          call rel (lsthed(i))
        end if
      end do
    end subroutine print_cross_reference ! of productions and states

  end subroutine PNTSET

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Print_Set

! $Log$
! Revision 1.6  2019/07/09 22:40:52  vsnyder
! Spiff some output
!
! Revision 1.5  2018/04/17 23:09:06  vsnyder
! J actually was being used as a global variable, but that's confusing, so
! pass it around as an argument.
!
! Revision 1.4  2018/04/17 22:40:54  vsnyder
! Declare local DO variables
!
! Revision 1.3  2014/01/14 00:11:42  vsnyder
! Revised LR completely
!
! Revision 1.2  2013/11/27 01:33:47  vsnyder
! Stop with stop-code 1 if not LR
!
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
