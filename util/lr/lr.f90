! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program LR

! Grammar analysis and parser table generation.

! The only known system dependent assumption is that all characters
! can be represented by integers in the range 0...255.  If this is
! not the case, change the declaration for TOGGLE in TOGCOM.

! A subprogram INIT, that does nothing in the portable version, is
! provided, in which one may perform system dependent initialization
! such as reading the command line and opening files.

! This is essentially a re-writing of the program LR, originally
! written by Alfred Shannon at the University of California Davis
! campus, under the direction of Charles Wetherell.  At the time of
! this writing, Shannon was at IBM Poughkeepsie.

! This version was written by:
!      William V. Snyder
!      Jet Propulsion Laboratory, Mail Stop 171-249
!      4800 Oak Grove Drive
!      Pasadena, CA 91109
!      818/354-6271, FTS/792-6271

! The input format for grammars is basically free-form.  In addition
! there is a punctuation convention that provides delimiters for
! productions and control for the program.  Almost any combination
! of characters may be used to construct vocabulary entries.

! The grammar should be translated to and input in a form similar to
! Backus-Naur form (BNF), the essential difference being in the
! punctuation.  The vocabulary of the language is inferred from the
! grammar:  Any symbol that occurs only on the right side is a
! terminal symbol and all others are nonterminal symbols.  If there
! is only one nonterminal symbol that appears only on the left side
! it is taken as the goal symbol.  If there is more than one the
! grammar is erroneous (it has useless productions).  If there are
! none the first nonterminal seen is taken as the goal symbol.  All
! of the productions that define a single nonterminal must occur
! together.

! All punctuation and control uses the ampersand character '&'.
! Whenever this character appears it signals a change into a mode
! where the following character defines the meaning or action.
! Normally, the two character string will then be discarded and
! input scanning resumed.  Some of the punctuation pairs are:

!     &A:  End of an alternate -- signals the end of one alternate
!          definition of a nonterminal.
!     &P:  End of productions -- signals the end of a group of
!          alternative definitions.
!     &G:  End of grammar -- terminate reading this grammar.

! The three pairs &A, &P, and &G delimit the definitions of
! nonterminals.  They are effective when they occur blank-delimited
! outside a vocabulary item.  If they occur within a vocabulary item
! they have the value of the null string, and any punctuation effect
! they would have had will be lost.

!     &S signals that the following token, which must be an integer,
!        is a syntax action index to be associated with the right
!        side currently being scanned.

!     &V signals that the next pair of tokens provide a value for a
!        vocabulary token.  The first item should be a terminal,
!        and the second must be an integer.

!     &B is an escaped blank -- create an item consisting of exactly
!        one blank.

!     && is an escaped & -- create an item consisting of exactly one
!        ampersand.

!     &< is an escaped < -- create an item consisting of exactly one
!        < character.  The < character otherwise has a special
!        meaning described below.

!     &> is an escaped > -- create an item consisting of exactly one
!        > character.  The > character otherwise has a special
!        meaning described below.

! The four previous pairs allow for the inclusion of characters that
! would normally have a special meaning.

!     &C signals the end of a line.  The rest of the line may be
!        used for comments.

! All other characters following an & invert control toggles.
! Control toggles are a set of on off switches (one for each
! character except those mentioned above) in the analyzer.  Every
! time &X (where X is any character other than the above) is seen
! in the input the state of the associated toggle is inverted.  The
! toggles are examined during execution of the analyzer to determine
! the selection of options.  Unless otherwise specified, toggles are
! initially off.  The currently defined toggles are:

!     &I list the input file.
!     &L output tables (initially on).
!     &M print the parsing machine (initially on).
!     &N list the grammar neatly and then quit.
!     &R list the grammar neatly (initially on).
!     &X include a vocabulary cross reference (initially on).
!     &2 print some debugging information about the configuration
!        analysis phase.
!     &3 print some more debugging information about the
!        configuration analysis phase.

! Vocabulary items are primarily blank delimited.  Blanks are
! discarded until a non-blank character appears.  If the non-blank
! character is an ampersand it is treated as described above.  If it
! is anything other than an ampersand or an unescaped <, everything
! up to the next blank is taken as one vocabulary item.  If the non-
! blank character is an unescaped <, everything up to and including
! the next unescaped > forms a vocabulary item.  The first item of
! a set of productions is taken to be the left hand side.  The
! following items up to the end of the productions make up the
! alternative right hand sides.  Thus the BNF
!     <EXPR> ::= <EXPR> + <TERM>
!              ! <TERM>
! would be input as
!     <EXPR> <EXPR> + <TERM> &A <TERM> &P

! The control and punctuation items must not appear within
! vocabulary items, and must be separated from some of them by
! blanks.  Also, the only vocabulary item that may begin with & is
! just & itself (although & may be buried inside longer items).  The
! effects of &A, &P and &G are cumulative.  That is, &P implies &A
! and &G implies &P.  For safety sake, it is easiest to delimit all
! items by blanks.

! An example grammar might be

!     S E &P
!     E E + T &A
!       T &P &C Ignore this part of the line.
!     T P ** T &A
!       P &P
!     P <Identifier> &A
!       ( E ) &G

! S is the goal symbol; S, E, T and P are nonterminals; +, **,
! <Identifier>, ( and ) are terminals.

! There will be a zero'th production

! <GOAL> ::= SOG S EOG

! Added for purposes of analysis.  This production is partially
! defined by block data.  The rest of the definition is created by
! the initialization phase of RDGRAM and in FNDGOL.  The symbol
! <GOAL> should not be used as a vocabulary item in the submitted
! grammar.

  use Analysis, only: Analyz
  use Chain_Context_Lists, only: CHNCSL
  use Connect, only: CONECT
  use Cross_Reference, only: XREF
  use Find_Goal, only: FNDGOL
  use Generate_Table, only: GENTAB
  use GRAMMAR, only: RDGRAM
  use Grounded, only: GROUND
  use LISTS, only: LISTS_INIT
  use Print_Grammar, only: PRNTGM
  use Print_Set, only: PNTSET
  use MD_INIT_M, only: MD_INIT
  use S3, only: SORTGM
  use TOGGLES, only: TOGGLE

  implicit NONE

!---------------------------- RCS Module Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  integer :: LNADQT = 1 ! Last inadequate state number, zero if grammar is LR.
                        ! Initially nonzero in case the grammar isn't analyzed.

  call md_init
  call lists_init
  call rdgram
  call fndgol
  call sortgm
  call conect
  call ground

  ! Toggle N (initially off) lists the grammar neatly and then quits.
  ! Toggle R (initially on) lists the grammar neatly.

  if (toggle(ichar('N')) /=0 .or. toggle(ichar('R')) == 0 ) call prntgm

  ! Toggle X (initially on) causes a cross reference to be printed.

  if (toggle(ichar('X')) == 0) call xref
  if (toggle(ichar('N')) == 0) then
     call analyz
     call pntset ( lnadqt )

     ! Toggle L (initially on) outputs the tables.

     if (toggle(ichar('L')) == 0) then
        call chncsl
        call gentab
     end if
  end if

  if ( lnadqt /= 0 ) stop 1

end program LR

! $Log$
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
