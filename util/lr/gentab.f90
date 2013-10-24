! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Generate_Table

  implicit NONE
  private
  public :: GENTAB

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine GENTAB

    use ANACOM, only: IFINAL, INDBAS, LENCSL, LISTCS, LSETS, NCSETS, &
                      NXTRED, NXTTRN
    use IO, only: OUNIT, OUTPUT, TBUNIT
    use LISTS, only: ITEM, NEXT
    use S1, only: LENGTH, MOVSTR
    use S3, only: ACTION, PRDIND, PRODCN, VALUE, VOCAB
    use S5, only: BASIS, RED, TRAN
    use TABCOM, only: NTERMS, NUMPRD, NVOC
    use TOGGLES, only: TOGGLE

    implicit NONE

  ! Generate the tables as Fortran named constants, or in unformatted form.

  ! Named constants or DATA statements if Fortran output is selected
  ! (toggle "F" is zero -- on) are

  ! IFINAL: The final state
  ! NUMPRD: The number of productions
  ! NSTATE: The number of states
  ! LSETS:  The total number of elements of lookahead sets
  ! NLOOKS: The number of states needing lookahead sets
  ! NCSETS: The number of lookahead sets
  ! NTRANS: The number of transitions
  ! NTERMS: The number of terminal symbols
  ! NVOC:   The total number of vocabulary symbols

  ! ENT(1:nstate): Indexed by state number.  Index of vocabulary symbol
  ! used to enter a state.

  ! FRED(0:nstate): Indexed by state number.  First (actually last)
  ! reductions.  Lookahead sets for state S are NSET(FRED(S-1)+1:FRED(S)).

  ! FTRN(0:nstate): Indexed by state number.  First (actually last)
  ! transitions.  Transitions for state S are in TRAN(FTRN(S-1)+1:FTRN(S)).

  ! TRAN(1:ntrans): New state indices, indexed by FTRN, q.v.

  ! NSET(1:nlooks): Lookahead set indices, indexed by FRED, q.v.

  ! PROD(1:nlooks): Production numbers to reduce if lookahead is in the set,
  ! indexed by FRED, q.v.

  ! LSET(0:ncsets): Last elements of lookahead sets, indexed by NSET, q.v.
  ! For lookahead set L, elements are LS(LSET(L-1)+1:LSET(L)).

  ! LS(1:lsets): Elements of lookahead sets.

  ! LENS(1:numprd): lengths of productions' right-hand sides, indexed by
  ! PROD.

  ! LHS(1:numprd): Productions' left-hand sides, indexed by PROD.

  ! ACT(1:numprd): Action to do if production is reduced, indexed by PROD.

  ! VOC(1:nvoc): Vocabulary symbol indices, from &V specs in the grammar. 
  ! Negative if no mapping has been provided.

  ! Numeric data, if Fortran output is not selected (toggle "F" is not zero)
  ! are output in 10I5 format.

  ! The first 13 items are:
  !    The number of states = NSTATE
  !    The final state = IFINAL
  !    The initial state = 1 (always)
  !    The size of the vocabulary = NVOC
  !    The number of terminal symbols = NTERMS
  !    The length of the longest symbol (computed here)
  !    The number of lookahead set choices = (NXTRED-1)/2
  !    The number of transitions = NXTTRN - 1
  !    The number of lookahead sets = NCSETS
  !    The number of productions = NUMPRD
  !    The length of the longest right hand side (computed here)
  !    The total number of lookahead sets = LSETS
  !    The total length of all vocabulary items

  ! Following these 13 items, the vocabulary items are output one per
  ! line in (I5,A) format.

  ! After the vocabulary items, arrays of numbers are output.  Each
  ! array starts on a new line:
  !    The entry symbol array (known as ENT in the parser).
  !    the first (actually last) lookahead array (known as FRED).
  !    The first (actually last) transition array (known as FTRN).
  !    The transition vector (known as TRAN).
  !    The pointer to lookahead sets array (known as NLOOKS).
  !    The reduction array (known as PROD).
  !    The lookahead sets.  Each lookahead set starts a new line.
  !       The first number is the number of terminals in the lookahead
  !       set.  The rest are token indices for the members of the set.
  !    The vector of right hand side lengths (known as LEN).
  !    The vector of left hand sides (known as LHS).
  !    The vector of syntax actions (known as ACTION).
  !    The vector of token values (known as VALUE).

  ! The following two arrays are not necessary for the parser, but
  ! may be useful for printing error messages.

  !    The right hand sides of the productions.  Each production
  !       starts on a new line.  The left hand sides and right hand
  !       side lengths were output previously.
  !    The configurations, without contexts.  Each state starts on a
  !       new line.  The first integer is the number of productions.
  !       Then for each production the production number and dot
  !       position are output.

  ! *****     External References     ********************************

  ! LENGTH  calculates the length of a symbol table item.
  ! MOVSTR  moves a symbol table item to LINE.
  ! OUTPUT  writes a line of text on the output file.

  !     *****     Local Variables     ************************************

  ! HLDIDX  is the index in the HOLD array.
  ! HOLD    is used to collect 10 integers for output.
  ! I       is a loop induction variable and subscript.
  ! IEND    is the upper limit for I in some loops.
  ! IPTR    is a list pointer.
  ! J       is a loop induction variable and subscript.
  ! K       is a subscript.
  ! L       is a subscript.
  ! LENPTR  is a pointer to an element of a list of lengths.
  ! LINE    is used for text assembly.
  ! LNGPRD  is the length of the longest right hand side of a production.
  ! LNGVCB  is the length of the longest vocabulary item.
  ! LSTPTR  is a pointer to a list element.
  ! TVOC    is the total length of all vocabulary items

    integer HLDIDX, HOLD(10), I, IEND, IPTR, J, K, L, LENPTR, N, NLOOKS, NSTATE
    character(len=120) :: LINE
    integer LNGPRD, LNGVCB, LSTPTR, TVOC
    integer :: LS(lsets), LSET(lsets+1:ncsets+lsets+1)
    !                                  nlooks       nstate          ntrans
    integer :: Work(0:max(lsets,ncsets,(nxtred-1)/2,indbas/5,nterms,nxttrn-1, &
                   &      numprd,nvoc))

  ! *****     Procedures     *****************************************

  ! Route further output to the table unit.

    ounit = tbunit

  ! output the first 13 numbers

    nlooks = (nxtred-1)/2         ! Number of lookahead set references
    nstate = indbas/5             ! Number of states
    if ( toggle(ichar('F')) /= 0 ) then ! Normally on
      hold(1) = nstate            ! Number of states
      hold(2) = ifinal            ! Final state
      hold(3) = 1                 ! Initial state
      hold(4) = nvoc              ! Number of vocabulary items
      hold(5) = nterms            ! Number of terminal symbols
      lngvcb = 0
      tvoc = 0
      do i = 1, nvoc
         lngvcb = max(lngvcb, length(vocab(i)))
         tvoc = tvoc + length(vocab(i))
      end do
      hold(6) = lngvcb            ! Longest vocabulary item
      hold(7) = nlooks            ! Number of lookahead set references
      hold(8) = nxttrn-1          ! Number of transitions
      hold(9) = ncsets            ! Number of lookahead sets
      hold(10) = numprd           ! Number of productions
      hldidx = 10
      call dump_the_hold_array
      lngprd = 0
      do i = 1, numprd
        lngprd = max(lngprd, prdind(i+1)-prdind(i)-1)
      end do
      hold(1) = lngprd            ! Longest right-hand-side
      hold(2) = lsets             ! Total size of all lookahead sets
      hold(3) = tvoc              ! Total size of the vocabulary
      hldidx = 3
      call dump_the_hold_array

    ! Output the vocabulary items, one per line, with the length in
    ! columns 1 - 5.

      do i = 1, nvoc
        l = 6
        write ( line(1:5), '(i5)' ) length(vocab(i))
        call movstr ( vocab(i), line, l, 80 )
        call output ( line(1:l-1) )
      end do
    else ! Output Fortran named constant declarations
      write ( line, 1 ) 'IFINAL', ifinal
    1 format ( 2x,"integer, parameter :: ", a, " = ", i0 )
      call output ( line )
      write ( line, 1 ) 'LSETS', lsets
      call output ( line )
      write ( line, 1 ) 'NCSETS', ncsets
      call output ( line )
      write ( line, 1 ) 'NLOOKS', nlooks
      call output ( line )
      write ( line, 1 ) 'NSTATE', nstate
      call output ( line )
      write ( line, 1 ) 'NTERMS', nterms
      call output ( line )
      write ( line, 1 ) 'NTRANS', nxttrn-1
      call output ( line )
      write ( line, 1 ) 'NUMPRD', numprd
      call output ( line )
      write ( line, 1 ) 'NVOC', nvoc
      call output ( line )
      hldidx = 0
    end if

  ! Output the entrance symbol array.  The entrance symbol for state
  ! 1 is output as 1 but state 1 is never entered so ENT(1) is never
  ! consulted.  (We use 1 instead of zero because zero is not a valid
  ! symbol number.  Some other processor may assume so.)

    n = 1
    work(1) = 1
    do i = 6, indbas-1, 5
      n = n + 1
      k = basis(i)
      work(n) = prodcn(prdind(basis(k))+basis(k+1)-1)
    end do
    call output_work ( Work(1:n), 'ENT', '(nstate)', 1, nstate )

  ! Output the first lookahead array (FRED).

    n = 0
    do i = 1, indbas+4, 5
      work(n) = (basis(i+4)+2)/2 - 1
      n = n + 1
    end do
    call output_work ( Work(0:indbas/5), 'FRED', '(0:nstate)', 0, nstate )

  ! Output the first transition array (FTRN)

    n = 0
    do i = 1, indbas+4, 5
      work(n) = basis(i+3) - 1
      n = n + 1
    end do
    call output_work ( Work(0:indbas/5), 'FTRN', '(0:nstate)', 0, nstate )

  ! Output the transition vector (TRAN).

    call output_work ( (tran(1: nxttrn - 1)+4)/5, &
                     & 'TRAN', '(ntrans)', 1, nxttrn - 1 )

  ! Output the pointers to lookahead vectors (NSET).

    call output_work ( item(red(2:nxtred-1:2)), &
                     & 'NSET', '(nlooks)', 1, (nxtred-1)/2 )

  ! Output the production index vector (PROD).

    call output_work ( red(1:nxtred-1:2), 'PROD', '(nlooks)', 1, (nxtred-1)/2 )

  ! Output the lookahead sets (LS and LSET).

    lenptr = lencsl
    lstptr = listcs
    if ( toggle(iachar('F')) == 0 ) then ! normally on
      i = lsets
      j = 0
      n = 0
    end if

    do while (lenptr > 0)
      iptr = next(item(lstptr))
      if ( toggle(iachar('F')) == 0 ) then ! normally on
        i = i + 1
        lset(i) = n
        n = n + item(lenptr)
        do while (iptr > 0)
          j = j + 1
          ls(j) = item(iptr)
          iptr = next(iptr)
        end do
        lset(i+1) = n
      else
        ! Each set starts a new line.
        ! The first number is the number of elements of the set.  Remaining
        ! elements are terminal token indices.
        hldidx = 1
        hold(1) = item(lenptr)
        do while (iptr > 0)
          if (hldidx >= 10) call dump_the_hold_array
          hldidx = hldidx + 1
          hold(hldidx) = item(iptr)
          iptr = next(iptr)
        end do
        call dump_the_hold_array
      end if
      lenptr = next(lenptr)
      lstptr = next(lstptr)
    end do

    if ( toggle(iachar('F')) == 0 ) then ! normally on
      call output_work ( lset, 'LSET', '(0:ncsets)', 0, ncsets )
      call output_work ( ls, 'LS', '(lsets)', 1, lsets )
    end if

  ! Output the right hand side length vector (LENS).

    call output_work ( prdind(2:numprd+1) - prdind(1:numprd) - 1, &
                     & 'LENS', '(numprd)', 1, numprd )

  ! Output the left hand side vector (LHS).

    call output_work ( prodcn(prdind(1:numprd)), 'LHS', '(numprd)', 1, numprd )

  ! Output the ACTION vector.

    call output_work ( action(1:numprd), 'ACT', '(numprd)', 1, numprd )

  ! Output the token value array (VALUE).

    call output_work ( value(1:nvoc), 'VAL', '(nvoc)', 1, nvoc )

  ! Output the right hand sides of the productions.

    do i = 1, numprd
      do l = prdind(i)+1, prdind(i+1)-1, 10
        do j = l, min(prdind(i+1)-1, l+9)
          hldidx = hldidx + 1
          hold(hldidx) = prodcn(j)
        end do
        call dump_the_hold_array
      end do
    end do

  ! Output the bases.

    do i = 1, indbas-1, 5
      hold(1) = (basis(i)-basis(i+5))/3
      hldidx = 1
      do l = basis(i), basis(i+5)+3, -3
        if (hldidx .ge. 10) call dump_the_hold_array
        hldidx = hldidx + 1
        hold(hldidx) = basis(l)
        if (hldidx .ge. 10) call dump_the_hold_array
        hldidx = hldidx + 1
        hold(hldidx) = basis(l+1)
      end do
      call dump_the_hold_array
    end do

    return

  contains

    subroutine dump_the_hold_array ( N, All )
      integer, intent(in), optional :: N, All
      integer :: I, J, K
      if ( toggle(iachar('F')) == 0 ) then ! normally on
        if ( present(n) ) then
          k = 3
          do i = 1, hldidx, 10
            do j = 1, min(i+9,hldidx)
              write ( line(k:k+4), '(i5)' ) hold(j)
              k = k + 5
              if ( j < min(i+9,hldidx) ) then
                line(k:k+1) = ', '
                k = k + 2
              end if
            end do
            if ( n /= all ) then
              line(k:k+3) = ', &'
            else
              line(k:k+2) = ' ]'
            end if
            call output ( line )
          end do
        end if
      else
        write ( line(1:5*hldidx), '(24i5)' ) hold(1:hldidx)
        call output (line(1:5*hldidx))
      end if
      hldidx = 0
    end subroutine dump_the_hold_array

    subroutine Gen_Parameter ( Name )
      character(len=*), intent(in) :: Name
      if ( toggle(iachar('F')) == 0 ) &
        & write ( line, '("  integer, parameter :: ",a,a)' ) trim(name), ' = [ &'
      call output ( line )
    end subroutine Gen_Parameter

    subroutine Output_Work ( Work, Name, Dim, First, Last )
      ! Dump Work(first:last) with Name, and Dim if as a named constant
      integer, intent(in) :: First, Last
      integer, intent(in) :: Work(First:)
      character(len=*), intent(in) :: Name, Dim
      integer, parameter :: MaxLines = 255 ! Fortran 2008 continuation limit
      integer :: I, J, K
      if ( toggle(iachar('F')) == 0 ) then ! normally on
        if ( last-first+1 <= 10*maxLines ) then ! Fits within one statement
          ! Output as a named constant
          line = '  integer, parameter :: ' // trim(name) // trim(dim) // ' = [ &'
          call output ( line )
          do i = first, last, 10
            k = 3
            do j = i, min(last,i+9)
              write ( line(k:k+4), '(i5)' ) work(j)
              k = k + 5
              if ( j < min(last,i+9) ) then
                line(k:k+1) = ', '
                k = k + 2
              end if
            end do
            if ( j > last ) then
              line(k:k+1) = ' ]'
              k = k + 2
            else
              line(k:k+2) = ', &'
              k = k + 3
            end if
            call output ( line(1:k) )
          end do
        else ! Too many continuation lines needed
          ! Output as a DATA statement
          line = '  integer :: ' // trim(name) // trim(dim)
          call output ( line )
          do i = first, last, 3
            write ( line, 1 ) (trim(name), j, work(j), j = i, min(i+2,last))
          1 format ( 2x,'data ', 3(a4, '(', i5, ') / ', i5, '/' : ', ') )
            call output ( line )
          end do
        end if
      else ! Output using (10i5)
        do i = first, last, 10
          write ( line(1:50), '(10i5)' ) work(i:min(i+9,last))
          call output ( line )
        end do
      end if
    end subroutine Output_Work

  end subroutine GENTAB

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Generate_Table

! $Log$
