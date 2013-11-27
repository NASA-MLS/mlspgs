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

    use ANACOM, only: IFINAL, INDBAS
    use Delete, only: DELCS
    use Error_Handler, only: Error
    use IO, only: OUTPUT
    use LISTS, only: ITEM, LINT, NEW, NEXT, REL
    use Print_List, only: PNTLST
    use S1, only: LENGTH, MOVSTR
    use S3, only: PRDIND, PRODCN, VOCAB
    use S5, only: BASIS, RED, TRAN
    use TABCOM, only: NUMPRD
    use TOGGLES
    implicit NONE

    integer, intent(out) :: LNADQT
    ! LNADQT  is the state number of the last inadequate state.  Zero if
    !         the grammar is adequate, it is LR(1).

    ! Print the configuration sets.

    ! *****     External References     ************************************

    ! DELCS   deletes a reference to a context set.
    ! ERROR   prints error messages.
    ! LENGTH  calculates the length of a vocabulary item.
    ! LINT    detects intersections between pairs of lists.
    ! MOVSTR  moves a vocabulary item from the symbol table.
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
    !         being detected, and subscript for LINE.
    ! LINK    is a pointer to a list node when cross reference lists are
    !         being constructed or printed.
    ! LINE    is used for message assembly.
    ! LSTHED  is an array of pointers to lists of states in which a
    !         a production is analyzed.
    ! REDUCE  is .TRUE. iff the last configuration printed was a reducing
    !         configuration, in which case the context is not printed (since
    !         it is the same as the lookahead).

    logical ADEQUT
    integer I, IBASE, IPR, IPTR, ITEMP, J, JDOT, JEND, JSTART
    integer K, L, LINK
    character(len=120) :: LINE
    integer LSTHED(NUMPRD)
    logical REDUCE

  !     *****     Procedures     *****************************************

    if (toggle(ichar('M')) /= 0) then
      call clean_up ! redundant context set references
      return
    end if
    lnadqt = 0
    call new (itemp)
    line = '1'
    line(39:81) = 'T H E   P A R S I N G   A U T O M A T O N'
    call output (line(1:81))
    line(1:14) = '0 State number'
    call output (line(1:14))
    line(1:26) = '     .  Production numbers'
    call output (line(1:26))
    line(1:27) = '     .       State Contents'
    call output (line(1:27))
    do i = 1, indbas-1, 5
      ipr = basis(i)
      call print_the_basis ! for the state
      call print_the_transitions ! from the state
      call print_the_reductions ! from the state
    end do
    call rel (itemp)
    if (lnadqt /= 0) then
      call error ('This grammar is not LR(1)',1)
      line = ' Last inadequate state:'
      write ( line(24:29), '(i6)' ) lnadqt
      call output (line(1:29))
    end if
    if (toggle(ichar('X')) == 0) &
      call print_cross_reference ! of productions and states

  contains

    subroutine clean_up ! redundant context set references
      do i = 1, indbas-1, 5
        ipr = basis(i)
        if (basis(ipr) == 1 .and. basis(ipr+1) > 3) ifinal = (i+4)/5
        do j = ipr, basis(i+5)+3, -3
          call delcs (basis(j+2))
        end do
      end do
    end subroutine clean_up ! redundant context set references

    subroutine print_the_basis ! for the state
      write ( line(1:6), '("0",i5)' ) (i+4)/5
      if (basis(ipr) == 1 .and. basis(ipr+1) > 3) ifinal = (i+4)/5
      do j = ipr, basis(i+5)+3, -3
        call print_the_production ! with a dot in the right side
        if (.not.reduce) then
          iptr = next(basis(j+2))
          if (iptr /= 0 .and. toggle(ichar('2')) /= 0) then
            line(14:21) = 'Context:'
            call pntlst (iptr, line, 23, 25)
          end if
        end if
        call delcs (basis(j+2))
      end do
    end subroutine print_the_basis ! for the state

    subroutine print_the_production ! with a dot in the right side

      ! Print the left side.

      write ( line(8:12), '(i5)' ) basis(j)
      l = 14
      ibase = prdind(basis(j))
      call movstr (vocab(prodcn(ibase)), line, l, 120)
      line(l+1:l+2) = '->'
      l = l + 4

      ! Print the right side of the production with a dot before the
      ! IDOT'th item.  IDOT is in BASIS(I+1).

      jdot = 1
      do k = ibase+1, prdind(basis(j)+1)-1
        if (length(vocab(prodcn(k)))+l > 120) then
          call output (line(1:l-1))
          l = 17
        end if
        if (jdot == basis(j+1)) then
          line(l:l) = '.'
          l = l + 2
        end if
        call movstr (vocab(prodcn(k)), line, l, 120)
        l = l + 1
        jdot = jdot + 1
      end do
      reduce = jdot == basis(j+1)
      if (reduce) then
        line(l:l) = '.'
        l = l + 2
      end if
      call output (line(1:l-1))
    end subroutine print_the_production ! with a dot in the right side

    subroutine print_the_transitions ! from the state
      jstart = basis(i+3)
      jend = basis(i+8) - 1
      if (jstart <= jend) then
        line(14:25) = 'Transitions:'
        l = 26
        do while (jstart <= jend)
          if (l > 115) then
            call output (line(1:l-1))
            l = 16
          end if
          write ( line(l:l+4), '(i5)' ) (tran(jstart)+4)/5
          l = l + 5
          jstart = jstart + 1
        end do
        call output (line(1:l-1))
      end if
    end subroutine print_the_transitions ! from the state

    subroutine print_the_reductions ! from the state
      adequt = .true.
      jstart = basis(i+4)
      jend = basis(i+9) - 2
      if (jstart <= jend) then
        line(14:24) = 'Reductions:'
        do j = jstart, jend, 2
          l = 26
          write ( line(l:l+4), '(i5)' ) red(j)
          l = l + 6
          iptr = next(red(j+1))
          call pntlst (iptr, line, l, 32)
          call check_for_conflicts ! in this state
        end do
        if (.not. adequt) then
          if (lnadqt /= 0) then
            line = ' Last inadequate state:'
            write ( line(24:29), '(i6)' ) lnadqt
            call output (line(1:29))
          end if
          lnadqt = (i+4)/5
        end if
      end if
    end subroutine print_the_reductions ! from the state

    subroutine check_for_conflicts ! in this state
      integer :: P ! position in LINE
      do l = basis(i+3), basis(i+8)-1
        k = basis(tran(l))
        item(itemp) = prodcn(prdind(basis(k))+basis(k+1)-1)
        if ( lint(itemp, next(red(j+1))) ) then
          adequt = .false.
          line = ' *** Intersection with transition to'
          write ( line(38:46), '(i0," on ")' ) (tran(l)+4)/5
          p = len_trim(line)+2
          call movstr ( vocab(item(itemp)), line, p, 80 )
          call output (line(1:p))
        end if
      end do
      if (j /= jstart) then
        do l = jstart, j-2, 2
          if ( lint(next(red(l+1)), next(red(j+1))) ) then
            adequt = .false.
            line = ' *** Intersection with reduction of'
            write ( line(36:40), '(i5)' ) red(l)
            call output (line(1:40))
          end if
        end do
      end if
    end subroutine check_for_conflicts ! in this state

    subroutine print_cross_reference ! of productions and states
      lsthed(1:numprd) = 0
      do i = 1, indbas-1, 5
        do j = basis(i), basis(i+5)+3, -3
          call new (link)
          item(link) = (i+4)/5
          next(link) = lsthed(basis(j))
          lsthed(basis(j)) = link
        end do
      end do

      ! Print the cross reference lists.

      line(1:1) = '1'
      line(36:84)='A U T O M A T O N   C R O S S   R E F E R E N C E'
      call output (line(1:84))
      line(1:11) = '0Production'
      call output (line(1:11))
      line(1:27) = '    .    Analyzed in states'
      call output (line(1:27))
      do i = 1, numprd
        link = lsthed(i)
        if (link /= 0) then
          write ( line(1:5), '(i5)' ) i
          l = 6
          do while (link /= 0)
            if (l > 115) then
            end if
            write ( line(l:l+4), '(i5)' ) item(link)
            l = l + 5
            link = next(link)
          end do
          call output (line(1:l-1))
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
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
