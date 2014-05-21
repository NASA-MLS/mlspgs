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

  subroutine GENTAB ( TBUNIT, Strings )

    use Basis_m, only: BASIS, IFINAL, INDBAS, Items
    use Chain_Context_Lists, only: LENCSL, LISTCS, LSETS, NCSETS                     
    use Declare_Vocabulary_m, only: Total
    use Output_m, only: Blanks, NewLine, Output, OutputOptions
    use LISTS, only: LIST
    use String_Table, only: Display_String, Get_String, String_Length
    use Tables, only: ACTION => Actions, First_Terminal, Last_Nonterminal, &
      & NTERMS, NUMPRD, PRDIND => Prod_Ind, PRODCN => Productions, VOCAB, &
      & Vocab_Names
    use Transitions_And_Reductions, only: NXTRED, NXTTRN, Red, Tran

    implicit NONE

    integer, intent(in) :: TBUNIT
    integer, intent(in) :: Strings(:) ! String indices of vocabulary names

  ! Generate the tables as Fortran named constants.

  ! IFINAL: The final state
  ! NUMPRD: The number of productions
  ! NSTATE: The number of states
  ! LSETS:  The total number of elements of lookahead sets
  ! NLOOKS: The number of states needing lookahead sets
  ! NCSETS: The number of lookahead sets
  ! NTRANS: The number of transitions
  ! NTERMS: The number of terminal symbols
  ! NVOC:   The total number of vocabulary symbols
  ! TOTAL:  The total number of symbols in all productions
  ! NBASPROD: Number of productions in all states of the grammar

  ! ENT(1:nstate): Indexed by state number.  Index of vocabulary symbol
  !    used to enter a state.

  ! FRED(0:nstate): Indexed by state number.  First (actually last)
  !    reductions.  Lookahead sets for state S are NSET(FRED(S-1)+1:FRED(S)).

  ! FTRN(0:nstate): Indexed by state number.  First (actually last)
  !    transitions.  Transitions for state S are in TRAN(FTRN(S-1)+1:FTRN(S)).

  ! TRAN(1:ntrans): New state indices, indexed by FTRN, q.v.

  ! NSET(1:nlooks): Lookahead set indices, indexed by FRED, q.v.

  ! PROD(1:nlooks): Production numbers to reduce if lookahead is in the set,
  !    indexed by FRED, q.v.

  ! LSET(0:ncsets): Last elements of lookahead sets, indexed by NSET, q.v.
  !    For lookahead set L, elements are LS(LSET(L-1)+1:LSET(L)).

  ! LS(1:lsets): Elements of lookahead sets.

  ! LENS(1:numprd): lengths of productions' right-hand sides, indexed by
  !    PROD.

  ! LHS(1:numprd): Productions' left-hand sides, indexed by PROD.

  ! ACT(1:numprd): Action to do if production is reduced, indexed by PROD.

  ! VAL(1:nvoc): Vocabulary type indices, from &V specs in the grammar. 
  !    Negative if no mapping has been provided.

  ! TEXT(1:nvoc): Vocabulary text indices.

  ! The following arrays are not necessary for the parser, but
  !    may be useful for printing error messages.

  ! PROD_IND(1:numprd+1): Production(i) is represented in
  !    PRODUCTIONS(prod_ind(i):prod_ind(i+1)-1).

  ! PRODUCTIONS(1:total): The productions, where each element is an index
  !    in the VAL and TEXT array.

  ! INDBAS(1:nstate+1): Indices of bases for each state.  Each element is
  !    the starting position in BASPRODS and DOTS of the production number
  !    and dot position for a basis.

  ! BASPRODS(1:nBaseProd): Production numbers in the bases

  ! DOTS(1:nBaseProd): Item in RHS of production before which a dot appears
  !    in the basis.

  ! *****     External References     ********************************

  ! OUTPUT  writes a line of text on the output file.

  !     *****     Local Variables     ************************************

  ! I       is a loop induction variable and subscript.
  ! IPTR    is a list pointer.
  ! J       is a loop induction variable and subscript.
  ! K       is a subscript.
  ! LENPTR  is a pointer to an element of a list of lengths.
  ! LINE    is used for text assembly.
  ! LINWID  is the maximum width of lines in Fortran output
  ! LSTPTR  is a pointer to a list element.

    integer :: I, IPTR, J, K, LENPTR, N, NLOOKS, NSTATE, &
      & NTRANS, NVOC
    character(len=120) :: LINE
    integer, parameter :: Linwid = 73
    integer :: LSTPTR
    integer :: LS(lsets), LSET(lsets+1:ncsets+lsets+1)
    integer :: NBaseProds ! Number of productions in the bases
    integer :: Old_Unit
    !                                  nlooks   nstate          ntrans
    integer :: Work(0:max(lsets,ncsets,nxtred-1,indbas-1,nterms,nxttrn-1, &
                   &      numprd,last_nonterminal - first_terminal + 1))

  ! *****     Procedures     *****************************************

  ! Route further output to the table unit.

    old_unit = outputOptions%prunit
    outputOptions%prunit = tbunit

    nBaseProds = basis(indbas)%item - 1
    call output ( '! This file is generated by the LR parser generator', advance='yes' )
    call output ( '! DO NOT EDIT THIS FILE MANUALLY', advance='yes' )
    call newLine
    call output ( '! Include this file in a module, after', advance='yes' )
    call output ( '! use Symbol_Types  ! to get token names', advance='yes' )
    call output ( '! use Tree_Types    ! to get tree node names', advance='yes' )
    call output ( '! implicit NONE', advance='yes' )
    call output ( '! public', advance='yes' )

    ! Output Fortran named constant declarations
    nlooks = nxtred - 1           ! Number of lookahead set references
    nstate = indbas - 1           ! Number of states
    ntrans = nxttrn-1
    nvoc = last_nonterminal - first_terminal + 1

  ! output the first 13 numbers

    write ( line, 1 ) 'IFINAL', ifinal
  1 format ( 2x,"integer, parameter :: ", a, " = ", i0 )
    call output ( line, advance='yes' )
    write ( line, 1 ) 'LSETS', lsets
    call output ( line, advance='yes' )
    write ( line, 1 ) 'NCSETS', ncsets
    call output ( line, advance='yes' )
    write ( line, 1 ) 'NLOOKS', nlooks
    call output ( line, advance='yes' )
    write ( line, 1 ) 'NSTATE', nstate
    call output ( line, advance='yes' )
    write ( line, 1 ) 'NTERMS', nterms
    call output ( line, advance='yes' )
    write ( line, 1 ) 'NTRANS', ntrans
    call output ( line, advance='yes' )
    write ( line, 1 ) 'NUMPRD', numprd
    call output ( line, advance='yes' )
    write ( line, 1 ) 'NVOC', nvoc
    call output ( line, advance='yes' )
    write ( line, 1 ) 'TOTAL', total
    call output ( line, advance='yes' )
    write ( line, 1 ) 'NBASPROD', nBaseProds
    call output ( line, advance='yes' )

  ! Output the entrance symbol array.  The entrance symbol for state
  ! 1 is output as 1 but state 1 is never entered so ENT(1) is never
  ! consulted.  (We use 1 instead of zero because zero is not a valid
  ! symbol number.  Some other processor may assume so.)

    work(1) = 1
    do i = 2, indbas-1
      k = basis(i)%item
      work(i) = prodcn(prdind(items(k)%prod)+items(k)%dot-1)
    end do
    call output_work ( Work(1:indbas-1), 'ENT', '(nstate)', 1, nstate )

  ! Output the first lookahead array (FRED).

    do i = 1, indbas
      work(i-1) = basis(i)%red - 1
    end do
    call output_work ( Work(0:indbas-1), 'FRED', '(0:nstate)', 0, nstate )

  ! Output the first transition array (FTRN)

    do i = 1, indbas
      work(i-1) = basis(i)%tran - 1
    end do
    call output_work ( Work(0:indbas-1), 'FTRN', '(0:nstate)', 0, nstate )

  ! Output the transition vector (TRAN).

    call output_work ( tran(1:ntrans), &
                     & 'TRAN', '(ntrans)', 1, nxttrn - 1 )

  ! Output the pointers to lookahead vectors (NSET).

    call output_work ( list(red(1:nlooks)%set)%item, &
                     & 'NSET', '(nlooks)', 1, nlooks )

  ! Output the production index vector (PROD).

    call output_work ( red(1:nlooks)%prod, 'PROD', '(nlooks)', 1, nlooks )

  ! Output the lookahead sets (LS and LSET).

    lenptr = lencsl
    lstptr = listcs
    i = lsets
    j = 0
    n = 0

    do while (lenptr > 0)
      iptr = list(list(lstptr)%item)%next
      i = i + 1
      lset(i) = n
      n = n + list(lenptr)%item
      do while (iptr > 0)
        j = j + 1
        ls(j) = list(iptr)%item
        iptr = list(iptr)%next
      end do
      lset(i+1) = n
      lenptr = list(lenptr)%next
      lstptr = list(lstptr)%next
    end do

    call output_work ( lset, 'LSET', '(0:ncsets)', 0, ncsets )
    call output_work ( ls, 'LS', '(lsets)', 1, lsets )

  ! Output the right hand side length vector (LENS).

    call output_work ( prdind(2:numprd+1) - prdind(1:numprd) - 1, &
                     & 'LENS', '(numprd)', 1, numprd )

  ! Output the left hand side vector (LHS).

    call output_work ( prodcn(prdind(1:numprd)), 'LHS', '(numprd)', 1, numprd )

  ! Output the action array ACT.

    call output_work ( action(1:numprd), 'ACT', '(numprd)', 1, numprd, which=2 )

  ! Output the token value array (VAL).

    call output_work ( vocab_names(1:nterms), 'VAL', '(nterms)', 1, nterms, &
      & which=3 )

  ! Output the indices where productions begin (PROD_IND).

    call output_work ( prdind(1:numprd+1), 'PROD_IND', '(numprd+1)', 1, numprd+1 )

  ! Output the productions (PRODUCTIONS).

    call output_work ( prodcn(1:total), 'PRODUCTIONS', '(total)', 1, total )

  ! Output the indices for the bases (INDBAS).

    call output_work ( basis(1:indbas)%item, 'INDBAS', '(nstate+1)', 1, nstate+1 )

  ! Output the production numbers for the bases (BASPRODS).

    call output_work ( items(basis(1)%item:basis(indbas)%item-1)%prod, &
      & 'BASPRODS', '(nbasprod)', 1, nBaseProds )

  ! Output the dot positions in the bases (DOTS).

    call output_work ( items(basis(1)%item:basis(indbas)%item-1)%dot, &
      & 'DOTS', '(nbasprod)', 1, nBaseProds )

  ! Generate a subroutine that fills the parser table.

    call output ( 'contains', advance='yes' )
    call output ( '  subroutine Init_Parser_Table ( T )', advance='yes' )
    call output ( '  ! Fill the parser table T.', advance='yes' )
    call output ( '    use Parser_Table_m, only: Parser_Table_t', advance='yes' )
    call output ( '    use Symbol_Table, only: Enter_Terminal', advance='yes' )
    call output ( '    type(parser_table_t), intent(out) :: T', advance='yes' )
    call output ( '    t%ifinal = ifinal', advance='yes' )
    call output ( '    t%lsets = lsets', advance='yes' )
    call output ( '    t%ncsets = ncsets', advance='yes' )
    call output ( '    t%nlooks = nlooks', advance='yes' )
    call output ( '    t%nstate = nstate', advance='yes' )
    call output ( '    t%nterms = nterms', advance='yes' )
    call output ( '    t%ntrans = ntrans', advance='yes' )
    call output ( '    t%numprd = numprd', advance='yes' )
    call output ( '    t%nvoc = nvoc', advance='yes' )
    call output ( '    t%total = total', advance='yes' )
    call output ( '    t%nBasProd = nBasProd', advance='yes' )
    call output ( '    allocate ( t%ent(1:nstate), source=ent )', advance='yes' )
    call output ( '    allocate ( t%fred(0:nstate), source=fred )', advance='yes' )
    call output ( '    allocate ( t%ftrn(0:nstate), source=ftrn )', advance='yes' )
    call output ( '    allocate ( t%tran(1:ntrans), source=tran )', advance='yes' )
    call output ( '    allocate ( t%nset(1:nlooks), source=nset )', advance='yes' )
    call output ( '    allocate ( t%prod(1:nlooks), source=prod )', advance='yes' )
    call output ( '    allocate ( t%lset(0:ncsets), source=lset )', advance='yes' )
    call output ( '    allocate ( t%ls(1:lsets), source=ls )', advance='yes' )
    call output ( '    allocate ( t%lens(1:numprd), source=lens )', advance='yes' )
    call output ( '    allocate ( t%lhs(1:numprd), source=lhs )', advance='yes' )
    call output ( '    allocate ( t%act(1:numprd), source=act )', advance='yes' )
    call output ( '    allocate ( t%val(1:nterms), source=val )', advance='yes' )
    call output ( '    allocate ( t%text(1:nvoc) )', advance='yes' )
    call output ( '    allocate ( t%prod_ind(1:numprd+1), source=prod_ind )', advance='yes' )
    call output ( '    allocate ( t%productions(1:total), source=productions )', advance='yes' )
    call output ( '    allocate ( t%indbas(1:nstate+1), source=indbas )', advance='yes' )
    call output ( '    allocate ( t%basProds(1:nBasProd), source=basProds )', advance='yes' )
    call output ( '    allocate ( t%dots(1:nbasprod), source=dots )', advance='yes' )
    do i = 1, nTerms
      call output ( i, before='    t%text(' )
      if ( vocab_names(i) < 0 ) then
        call display_string ( strings(-vocab_names(i)), &
          & before=') = enter_terminal("_' )
      else
        call display_string ( vocab(i), before=') = enter_terminal("' )
      end if
      call output ( '",0)', advance='yes' )
    end do
    do i = nTerms+1, nvoc
      call output ( i, before='    t%text(' )
      call display_string ( vocab(i), before=') = enter_terminal("_' )
      call output ( '",0)', advance='yes' )
    end do
    call output ( '  end subroutine Init_Parser_Table', advance='yes' )

    outputOptions%prunit = old_unit ! Back to whichever unit OUTPUT was using

  contains

    subroutine Gen_Parameter ( Name )
      character(len=*), intent(in) :: Name
      write ( line, '("  integer, parameter :: ",a,a)' ) trim(name), ' = [ &'
      call output ( line )
    end subroutine Gen_Parameter

    subroutine Output_Numeric_Named_Value ( Work, Name, Dim, First, Last )
      ! Dump Work(first:last) with Name and Dim, as a named constant
      integer, intent(in) :: First, Last
      integer, intent(in) :: Work(First:)
      character(len=*), intent(in) :: Name, Dim
      integer :: I, J
      line = '  integer, parameter :: ' // trim(name) // trim(dim) // ' = [ &'
      call output ( line, advance='yes' )
      do i = first, last, 10
        call blanks ( 4 )
        do j = i, min(last,i+9)
          call output ( work(j), 5 )
          if ( j /= last ) call output ( ', ' )
        end do
        if ( j > last ) then
          call output ( ' ]', advance='yes' )
        else
          call output (  ' &', advance='yes' )
        end if
      end do
    end subroutine Output_Numeric_Named_Value

    subroutine Output_String ( Work, Name, Dim, First, Last, Width, Which )
      integer, intent(in) :: First, Last
      integer, intent(in) :: Work(First:)
      character(len=*), intent(in) :: Name, Dim
      integer, intent(in) :: Width ! of widest field
      integer, intent(in) :: Which ! 2 = ACT, 3 = VAL
      integer :: I, J
      character(10) :: Format
      logical :: More  ! Outputting ACT and Work(i) is not numeric
      integer :: NPL   ! Number Per Line
      character(width+3) :: Temp
      integer :: W     ! Width of one string

      call output ( &
        & '  integer, parameter :: ' // trim(name) // trim(dim) // ' = [ &', &
        & advance='yes' )
      npl = max(1,linwid/width)
      write ( format, '("(i",i0,")")' ) width-2 ! in case of numeric output
      do j = first, last, npl
        call blanks ( 4 )
        do i = j, min(j+npl-1,last)
          if ( work(i) == 0 ) then
            write ( line, format ) merge(0,-i,which==2)
            call output ( line(:width-2) )
          else if ( which == 2 ) then ! ACT
            w = string_length(abs(work(i)))
            call get_string ( abs(work(i)), temp )
            more = verify(temp(:w),'0123456789') /= 0 ! Not numeric
            if ( more ) then
              call blanks ( width - w - 7 ) ! -7 because width includes ", "
                                            ! "10*" before, and "+1" after
              call output ( '10*' )
              temp(w+1:w+2) = '+' // merge('1','2',work(i)>0)
              w = w + 2
            else
              call blanks ( width - w - 2 )
            end if
            call output ( temp(1:w) )
          else ! which == 3 => VAL
            if ( work(i) > 0 ) then
              w = string_length(work(i))
              call blanks ( width - w - 2 ) ! -2 because width includes ", "
              call display_string ( work(i) )
            else
              write ( line, format ) work(i)
              call output ( line(:width-2) )
            end if
          end if
          if ( i /= last ) call output ( ', ' )
        end do
        if ( i > last ) then
          call output ( ' ]', advance='yes' )
        else
          call output ( '&', advance='yes' )
        end if
      end do
    end subroutine Output_String

    subroutine Output_Work ( Work, Name, Dim, First, Last, Which )
      ! Dump Work(first:last) with Name, and Dim if as a named constant
      integer, intent(in) :: First, Last
      integer, intent(in) :: Work(First:)
      character(len=*), intent(in) :: Name, Dim
      integer, intent(in), optional :: Which
                                   ! If present and == 2, output symbol from
                                   ! symbol table, not integer from WORK.
                                   ! If present and == 3, output symbol from
                                   ! symbol table, preceded by "10*" and
                                   ! followed by "+1" if the string index is
                                   ! positive and "+2" if it's negative.
      integer, parameter :: MaxLines = 255 ! Fortran 2008 continuation limit
      integer :: I, J, K, W
      integer :: Width             ! of longest symbol
      integer :: My_First, My_Last
      character(31) :: My_Name, My_Dim
      w = 1
      if ( present(which) ) w = which
      width = 5   ! for numeric values, except for ACT and VAL
      if ( w /= 1 ) then
        ! Work out length of longest symbol to determine how many
        ! to put on each line
        width = 5
        do i = first, last
          if ( w == 2 ) then ! ACTIONS
            if ( work(i) /= 0 ) width = max(width,string_length(abs(work(i))),2)
          else ! Vocab_Names
            if ( work(i) > 0 ) width = max(width,string_length(work(i)),2)
          end if
        end do
        width = width+merge(7,2,which==2) ! Room for 10*name+? if ACT
      end if
      k = ( last-first ) / ( linwid/width * maxLines )
      if ( k == 0 ) then ! Fits within one statement
        ! Output as one named constant
        if ( w == 1 ) then
          call output_numeric_named_value ( work(first:), name, dim, first, &
            & last )
        else
          call output_string ( work(first:), name, dim, first, last, width, w )
        end if
      else
        ! Output as several named constants
        my_first = first
        do j = 0, k
          write ( my_name, '(a,"_",i0)' ) trim(name), j
          my_last = min(last, my_first + 10*maxLines - 1)
          write ( my_dim, '("(",i0,":",i0,")")' ) my_first, my_last
          if ( w == 1 ) then
            call output_numeric_named_value ( work(my_first:), my_name, &
              & my_dim, my_first, my_last )
          else
            call output_string ( work(my_first:), my_name, my_dim, &
              & my_first, my_last, width, w )
          end if
          my_first = my_last + 1
        end do
        line = '  integer, parameter :: ' // trim(name) // trim(dim) // ' = [ &'
        call output ( line )
        i = 5
        do j = 0, k
          write ( my_name, '(a,"_",i0)' ) trim(name), j
          if ( i + len_trim(my_name) > 78 ) then
            line(i:) = ', &'
            call output ( line )
            i = 5
          end if
          if ( line(i-1:i-1) /= '' ) then
            line(i:) = ', '
            i = i + 2
          end if
          line(i:) = my_name
          i = len_trim(line) + 1
        end do
        line(i+1:) = ' ]'
        call output ( line )
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
! Revision 1.4  2014/01/18 02:44:02  vsnyder
! NXTRED should have been NLOOKS in two places
!
! Revision 1.3  2014/01/14 00:11:42  vsnyder
! Revised LR completely
!
! Revision 1.2  2013/12/12 01:53:00  vsnyder
! Remove unused cruft
!
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
