! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Analysis

  implicit NONE
  private
  public :: ANALYZ

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine ANALYZ

    use Basis_m, only: ADDBAS, DEQUE, NEWBAS
    use Complete, only: Closures_Index, Closures, Complt
    use FIRST_SETS, only: FIND_FIRST      ! Find FIRST sets
    use LISTS, only: LIST, NEW
    use New_Context_Set, only: NEWCS
    use NULLABLE_M, only: FIND_NULLABLE   ! Find nullable nonterminals
    use Sort_Configurations, only: SORTCG
    use Tables, only: HEADCS, HEADEN, NTERMS, NVOC, PRDIND => Prod_Ind, &
      & PRODCN => Productions
    use Toggles, only: GEN, Toggle
    use Trace, only: Trace_Begin, Trace_End
    use Transitions_And_Reductions, only: TRNRED

    implicit NONE

    ! Analyze the grammar.  First initialize the BASIS array by putting
    ! a state corresponding to the first production with the dot before
    ! the first symbol into the BASIS array.  Then analyze each state
    ! fetched from a queue by completing the state (calculating its
    ! closure) and then constructing new states as necessary.

    ! *****     External References     ********************************

    ! ADDBAS  adds a new configuration to the basis.
    ! COMPLT  completes the basis by adding the closure.
    ! DEQUE   deques a basis for processing.
    ! NEWBAS  gets ready for construction of a new basis.
    ! NEWCS   constructs a new context set.
    ! SORTCG  Sorts a configuration.
    ! TRNRED  calculates transitions and reductions for the basis, and
    !         new states as necessary.

    ! *****     Local Variables     ************************************

    ! I       is a loop induction variable and subscript, and pointer to the
    !         basis currently being processed.
    ! J1, J2  are the bounds for the part of Closures to use
    ! N       is the pointer to the context set for production 1.
    ! NPTR    is the pointer to the initial context list for production 1.

    integer I, J1, J2, N, NPTR

    ! *****     Procedures     *****************************************

    if ( toggle(gen) ) call trace_begin ( 'Analysis' )
    call find_nullable
    call find_first

    allocate ( headcs(1:nterms), source = 0 )
    allocate ( headen(1:nvoc), source = 0 )

    call newbas ( i )
    call new ( nptr ) ! Initial context list for production 1.
    list(nptr)%item = prodcn(prdind(1) + 1) ! <SOG>
    call newcs ( nptr, n )

    !           iptr   ndot
    !               npr   nset
    call addbas ( i, 1, 1, n )         ! <GOAL> ::= . <SOG> start <EOG>

    do while ( i > 0 )
       call complt ( i )
       j1 = closures_index(i-1) + 1
       j2 = closures_index(i)
       call sortcg ( closures(j1:j2) )
       call trnred ( i, closures(j1:j2) )
       call deque ( i )
    end do
    if ( toggle(gen) ) call trace_end ( 'Analysis' )

  end subroutine ANALYZ

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Analysis

! $Log$
! Revision 1.3  2014/01/14 00:11:42  vsnyder
! Revised LR completely
!
! Revision 1.2  2013/12/12 01:53:00  vsnyder
! Remove unused cruft
!
! Revision 1.1  2013/10/24 22:41:13  vsnyder
! Initial commit
!
