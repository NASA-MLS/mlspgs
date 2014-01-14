! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Tables

  ! Tables used throughout the LR parser generator.

  implicit NONE
  public
  save ! The default since 2003

  integer :: First_Nonterminal ! Position in sorted list
  integer :: First_Terminal    ! Position in sorted list
  integer :: Last_Nonterminal  ! Position in sorted list

  integer :: NTERMS            ! Number of terminals
  integer :: NUMPRD            ! Number of productions
  integer :: NVOC              ! Number of symbols in the vocabulary

  integer, allocatable :: Actions(:)  ! String index, indexed by production
                               ! number, negative if "?" action

  integer, allocatable :: First_Production(:) ! First element of Productions
                               ! for a nonterminal; indexed by string index

  integer, allocatable :: HEADCS(:) ! Head of list of context sets for a
                               ! nonterminal symbol

  integer, allocatable :: HEADEN(:) ! Head of list of states entered by a
                               ! terminal symbol

  integer, allocatable :: Productions(:) ! LHS, RHS... indices in the
                               ! permutation table of sorted strings

  integer, allocatable :: Prod_Ind(:) ! Index in Productions of beginning of
                               ! each production; indexed by production number,
                               ! 1 + NUMPRD elements.

  integer, allocatable :: Vocab(:)    ! String indices of sorted vocabulary,
                               ! in order first by Terminal, Nonterminal,
                               ! Vocabulary name, and Action, then by string

  integer, allocatable :: Vocab_Names(:) ! String indices of terminal symbols
                               ! names, from "A = B" lines in the grammar, in
                               ! the same order as the terminals are in VOCAB,
                               ! if any, else the negative of the element's
                               ! subscript.

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

end module Tables

! $Log$
