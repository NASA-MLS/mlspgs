! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Print_The_Vocabulary_m

  ! Sort the vocabulary, first on symbol type (Terminal, Nonterminal,
  ! Vocabulary name, Action) and then on symbol text.
  ! Print the vocabulary in two columns, terminals on the left,
  ! nonterminals on the right.

  implicit NONE
  private

  public :: Print_The_Vocabulary

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Procedures     =====================================

  subroutine Print_The_Vocabulary ( P )
    ! Sort the vocabulary, first on type ( (Terminal, Nonterminal,
    ! Vocabulary name, Action) and then on symbol text.
    ! Then print the sorted vocabulary

    use Declaration_Table, only: Action, Declaration, Decls, Empty, Get_Decl, &
      & Nonterminal, Null_Decl, Terminal, Vocabulary
    use Output_m, only: Newline, Output
    use Processor_Dependent, only: NewPage
    use String_Table, only: Get_String, How_Many_Strings, String_Length
    use Tables, only: First_Nonterminal, First_Terminal, Last_Nonterminal, &
      & NTerms, NVoc
    use Toggles, only: Switches

    integer, allocatable, intent(out) :: P(:) ! Permutation vector for sorting the vocabulary

    type(decls) :: Decl
    integer :: First_Action, First_Undeclared, First_Vocab
    integer :: I, J
    character(120) :: Line
    integer :: N_Nonterminal, N_Terminal
    integer :: NNTerms          ! Number of nonterminals
    integer :: W                ! Width of a string
    logical :: Watch

    watch = index(switches,'dvoc') /= 0

    ! Sort the vocabulary
    allocate ( p(how_many_strings()) )
    call gsortp ( compar, size(p), p )

    ! Determine where the vocabulary type boundaries are
    first_action = huge(0)
    first_nonterminal = huge(0)
    first_terminal = huge(0)
    first_undeclared = huge(0)
    first_vocab = huge(0)
    do i = 1, size(p)
      decl = get_decl(p(i),nonterminal)
      if ( decl%type /= nonterminal ) decl = declaration(p(i))
      select case ( decl%type )
      case ( action )
        first_action = min(first_action,i)
      case ( nonterminal )
        first_nonterminal = min(first_nonterminal,i)
      case ( terminal )
        first_terminal = min(first_terminal,i)
      case ( vocabulary )
        first_vocab = min(first_vocab,i)
      case ( null_decl, empty )
        first_undeclared = min(first_undeclared,i)
      end select
    end do
    last_nonterminal = min(first_action-1,first_vocab-1,first_undeclared-1,size(p))

    n_nonterminal = min(first_vocab,first_undeclared) - first_nonterminal
    n_terminal = first_nonterminal - first_terminal

    nnterms = last_nonterminal - first_nonterminal + 1
    nterms = first_nonterminal - first_terminal
    nvoc = nnterms + nterms

    if ( watch ) then
      call output ( 'After revision by Print_The_Vocabulary:', advance='yes' )
      call output ( nnTerms, 5 ); call output ( ' Nonterminals', advance='yes' )
      call output ( nTerms, 5 );  call output ( ' Terminals', advance='yes' )
      call output ( nVoc, 5 );    call output ( ' Symbols', advance='yes' )
    end if

    ! Print terminals and nonterminals
    call output ( newPage, dont_asciify=.true. )
    line(1:50) = '     T E R M I N A L S'
    line(51:) = 'N O N T E R M I N A L S'
    call output ( trim(line), advance='yes' )
    call newline
    do i = 1, max(n_nonterminal,n_terminal)
      line = ''
      if ( i <= n_terminal ) then
        write ( line(1:4), '(i4)' ) i
        j = i + first_terminal - 1
        decl = get_decl(p(j),terminal)
        call get_string ( p(j), line(6:) )
        if ( decl%value /= 0 ) then
          w = string_length(p(j))
          line(6+w:9+w) = ' = '
          call get_string ( decl%value,line(9+w:) )
        end if
      end if
      if ( i <= n_nonterminal ) then
        j = i + first_nonterminal - 1
        write ( line(46:49), '(i4)' ) i + n_terminal
        call get_string ( p(j), line(51:) )
      end if
      call output ( trim(line), advance='yes' )
    end do

  end subroutine Print_The_Vocabulary

  integer function COMPAR ( I, J )

    use Declaration_Table, only: Declaration, Decls, Get_Decl, Nonterminal, &
      & Null_Decl
    use String_Table, only: Compare_Strings

    integer, intent(in) :: I, J ! String indices
    type(decls) :: I_Decl, J_Decl
    integer :: I_Type, J_Type   ! Type indices

    ! Get the string declarations, with priority given to nonterminals
    ! because they might at first have been declared as terminals
    i_decl = get_decl(i,nonterminal)
    if ( i_decl%type /= nonterminal ) i_decl = declaration(i)
    i_type = merge(i_decl%type,99,i_decl%type/=null_decl)
    j_decl = get_decl(j,nonterminal)
    if ( j_decl%type /= nonterminal ) j_decl = declaration(j)
    j_type = merge(j_decl%type,99,j_decl%type/=null_decl)
    compar = i_type - j_type
    if ( compar == 0 ) compar = compare_strings ( i, j, caseless=.false. )

  end function COMPAR

      subroutine GSORTP (COMPAR, N, P)
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  All rights reserved.  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!>> 1998-01-20  GSORTP  Snyder  Allow not initializing P.
!>> 1996-05-01  GSORTP  Krogh   Changes to use .C. and C%%.
!>> 1995-11-17  GSORTP  Krogh   Converted SFTRAN to Fortran 77.
!>> 1991-04-02  GSORTP  Snyder  Repair no permutation vector if m-n < 10
!>> 1988-11-22  GSORTP  Snyder  Initial code.
!
!     Sort an N-vector of objects of unknown type and organization.
!     P is set so that the P(J)'th element of the original sequence is
!     the J'th element of the sorted sequence.  The order is defined by
!     the integer function COMPAR.  An invocation COMPAR(I,J) should
!     return -1 if the I'th element of the data is to preceed the J'th
!     element in the sorted sequence, +1 if the J'th element is to
!     preceed the I'th element in the sorted sequence, and 0 if the I'th
!     and J'th elements are the same.
!
!     This subprogram is unaware of the data, and cannot manipulate it.
!     It is the caller's responsibility to make the data known to the
!     COMPAR function.
!
      interface
        integer function COMPAR ( I, J )
          integer, intent(in) :: I, J
        end function COMPAR
      end interface
      integer, intent(in) :: N
      integer, intent(out) :: P(*)
!
!     *****     Local Variables     ************************************
!
! BL      Left bound of the sub-array to be sorted at the next step.
! BR      Right bound of the sub array to be sorted at the next step.
! CL      Current left bound of the unsorted sub-array.
! CR      Current right bound of the unsorted sub-array.
! MYN     My N.
! PARTN   is the partition element.
! PTEMP   holds elements of P during exchanges.
! STACKL  keeps track of the left bounds of sub-arrays waiting to be
!         sorted.
! STACKR  keeps track of the right bounds of sub-arrays waiting to be
!         sorted.
! STKTOP  keeps track of the top of the stacks.
!
      integer BL,BR,CL,CR,MYN,PARTN,PTEMP
      integer STACKL(32),STACKR(32),STKTOP
!
!     *****     Executable Statements     ******************************
!
      do cl = 1, n
         p(cl)=cl
      end do
      myn = abs(n)
      if (myn >= 10) then
         stktop=1
         stackl(1)=1
         stackr(1)=myn
!           Start until loop
         do while ( stktop /= 0 )
            bl=stackl(stktop)
            br=stackr(stktop)
            stktop=stktop-1
!           Choose a partitioning element.  Use the median of the first,
!           middle and last elements.  Sort them so the extreme elements
!           can serve as sentinels during partitioning.
            cl=(bl+br)/2
            partn=p(cl)
            if (compar(p(bl),partn) > 0) then
               p(cl)=p(bl)
               p(bl)=partn
               partn=p(cl)
            end if
            if (compar(p(bl),p(br)) > 0) then
               ptemp=p(bl)
               p(bl)=p(br)
               p(br)=ptemp
            end if
            if (compar(partn,p(br)) > 0) then
               p(cl)=p(br)
               p(br)=partn
               partn=p(cl)
            end if
            p(cl)=p(br-1)
            p(br-1)=partn
!           Partition the sub-array around PARTN.  Exclude the above
!           considered elements from partitioning because they're al-
!           ready in the correct subfiles.  Stop scanning on equality to
!           prevent files containing equal values from causing a loop.
            cl=bl
            cr=br-1
            do
              do
                cl=cl+1
                if (compar(p(cl),partn) >= 0) exit
              end do
              do
                cr=cr-1
                if (compar(p(cr),partn) <= 0) exit
              end do
              if (cl > cr) exit
              ptemp=p(cl)
              p(cl)=p(cr)
              p(cr)=ptemp
            end do
!           Put sub-arrays on the stack if they're big enough.  Put the
!           larger under the smaller, so the smaller will be done next.
!           This makes the upper bound of the stack depth log2 (myn).
!           (The "Hibbard" modification of quicksort).
            if (cl-bl > br-cr) then
               if (cl-bl > 10) then
                  stktop=stktop+1
                  stackl(stktop)=bl
                  stackr(stktop)=cr
               end if
               if (br-cr > 10) then
                  stktop=stktop+1
                  stackl(stktop)=cl
                  stackr(stktop)=br
               end if
            else
               if (br-cr > 10) then
                  stktop=stktop+1
                  stackl(stktop)=cl
                  stackr(stktop)=br
               end if
               if (cl-bl > 10) then
                  stktop=stktop+1
                  stackl(stktop)=bl
                  stackr(stktop)=cr
               end if
            end if
!           End until loop
         end do ! while ( stktop /= 0 )
      end if
!     Clean up small subfiles by using insertion sort on everything.
      do cr = 2, myn
         ptemp=p(cr)
         cl=cr
         do while (compar(p(cl-1),ptemp) > 0)
            p(cl)=p(cl-1)
            cl=cl-1
            if (cl <= 1) exit
         end do
         p(cl)=ptemp
      end do

      end subroutine GSORTP

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Print_The_Vocabulary_m

! $Log$
! Revision 1.2  2019/07/09 20:27:24  vsnyder
! Compute number of nonterminals correctly if there is no vocabulary
!
! Revision 1.1  2014/01/14 00:15:05  vsnyder
! Initial commit of new module for new LR
!
