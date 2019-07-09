! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Declaration_Table

  use Tree, only: Null_Tree

  implicit NONE
  private

  ! Public procedures
  public :: Allocate_Decl, Deallocate_Decl, Declaration, Declare, &
            Declared, Dump_Decl, Dump_A_Decl, Dump_1_Decl, Get_Decl, &
            Init_Decl, Redeclare

  ! Public type
  public :: Decls

  ! Public declaration types
  public :: Action, Empty, Nonterminal, Null_Decl, Terminal, Type_Names, &
          & Vocabulary

  ! The type indices for declarations are such that when the symbol table
  ! is sorted, first on type and then on string value, the order will be
  ! NULL_DECL, Terminals, Nonterminals, Vocabulary, Actions, Empty

  integer, parameter :: NULL_DECL = 0   ! Index of the null declaration

  integer, parameter :: Terminal = 1    ! Symbol appears as a terminal
  integer, parameter :: Nonterminal = 2 ! Symbol appears as a nonterminal
  integer, parameter :: Vocabulary = 3  ! Symbol appears as a vocabulary value
  integer, parameter :: Action = 4      ! Symbol appears as an action
  integer, parameter :: Empty = 5       ! Not one of the previous

  integer, parameter :: Last_Type = Empty
  ! maxval(action, nonterminal, terminal, vocabulary )

  type :: Decls
    integer :: Type = empty
    integer :: Value = 0         ! String index of name of symbol in
                                 ! vocabulary definition
    integer :: Tree = null_tree  ! Tree node index related to declaration
    integer :: Prior = null_decl ! Prior declaration of same symbol
  end type Decls

  ! Strings need to be in the same order as type names
  character(len=*), parameter :: Type_Names(null_decl:last_type) = &
    [ 'Null_Decl  ', 'Terminal   ', 'Nonterminal', &
      'Vocabulary ', 'Action     ', 'Empty      ' ]

! -----     Private declarations     -----------------------------------
  type(decls), save, allocatable :: Decl_Table(:)
  integer, save :: Num_Decls = 0        ! amount of Decl_Table used
  integer, save, allocatable :: Symbol_Decl(:) ! indexed by string
                                        ! index, gives index in decl_table.

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains ! ====     Procedures     =====================================

! ------------------------------------------------  ALLOCATE_DECL  -----
  subroutine ALLOCATE_DECL ( NDECLS, STAT )
  ! Allocate NDECLS declarations.  Allocate STRING_TABLE_SIZE symbol
  ! declarations (indexes from STRING_TABLE to DECL_TABLE).  Also does
  ! INIT_DECL.
    use, intrinsic :: ISO_Fortran_ENV, only: Error_Unit
    use String_Table, only: String_Table_Size
    integer, intent(in) :: NDECLS       ! Number to allocate
    integer, intent(out), optional :: STAT   ! From ALLOCATE statement
    integer :: MY_STAT
    character(120) :: ERMSG ! From allocate
    if ( allocated(decl_table) ) deallocate(decl_table)
    if ( allocated(symbol_decl) ) deallocate(symbol_decl)
    allocate ( decl_table(0:ndecls), stat=my_stat, errmsg=ermsg )
    if ( my_stat /= 0 ) then
      if ( present(stat) ) then
        stat = my_stat
        return
      end if
      write ( error_unit, '(a,i0/a)' ) &
        & 'DECL_TABLE%ALLOCATE_DECL-E- Unable to allocate DECL_TABLE; STAT = ', &
        & my_stat, trim(ermsg)
      stop
    end if
    allocate ( symbol_decl(0:string_table_size()), stat=my_stat, errmsg=ermsg )
    if ( my_stat /= 0 ) then
      if ( present(stat) ) then
        stat = my_stat
        return
      end if
      write ( error_unit, '(a,i0/a)' ) &
        & 'DECL_TABLE%ALLOCATE_DECL-E- Unable to allocate SYMBOL_DECL; STAT = ', &
        & my_stat, trim(ermsg)
      stop
    end if
    call init_decl
    symbol_decl = null_decl

  end subroutine ALLOCATE_DECL

! ----------------------------------------------  DEALLOCATE_DECL  -----
  subroutine DEALLOCATE_DECL
    if ( allocated(decl_table) ) deallocate ( decl_table )
    if ( allocated(symbol_decl) ) deallocate( symbol_decl )
  end subroutine DEALLOCATE_DECL

! --------------------------------------------------  DECLARATION  -----
  type(decls) function DECLARATION ( STRING )
    ! Return the most recent declaration for STRING
    integer, intent(in) :: STRING  ! String index for which declaration needed
    if ( string > ubound(symbol_decl,1) ) & ! Assume string_table is increased
      & call increase_symbol_decl
    declaration = decl_table(symbol_decl(string))
  end function DECLARATION

! ------------------------------------------------------  DECLARE  -----
  subroutine DECLARE ( String, type, Tree, Value )
    use, intrinsic :: ISO_Fortran_ENV, only: Error_Unit
    use String_Table, only: Display_String, How_Many_Strings
    use Toggles, only: Levels, Tab, Toggle
    integer, intent(in) :: String  ! String index of name to declare
    integer, intent(in) :: Type    ! Type of object, e.g. UNITS, LABEL...
    integer, intent(in) :: Tree    ! Index of tree node of declaration
    integer, intent(in), optional :: Value

    integer :: MyValue, N
    integer :: STAT
    type(decls), allocatable :: New_Decl(:)
    character(120) :: ERMSG ! From allocate

    myValue = 0
    if ( present(value) ) myValue = value
    if ( .not. allocated(decl_table) ) call allocate_decl ( how_many_strings() )
    n = ubound(decl_table,1)
    num_decls = num_decls + 1
    if ( num_decls > n ) then
    ! Double size of declaration table
      allocate( new_decl(0:2*n), stat=stat, errmsg=ermsg )
      if ( stat /= 0 ) then
        write ( error_unit, '(a,i0/a)' ) &
          & 'DECL_TABLE%ALLOCATE_DECL-E- Unable to allocate NEW_DECL; STAT = ', &
          & stat, trim(ermsg)
        stop
      end if
      new_decl(0:n) = decl_table
      call move_alloc ( new_decl, decl_table )
    end if
    if ( string > ubound(symbol_decl,1) ) & ! Assume string_table is increased
      & call increase_symbol_decl

                              !     type  value    tree  prior
    decl_table(num_decls) = decls ( type, myValue, tree, symbol_decl(string) )
    symbol_decl(string) = num_decls
    if ( toggle(tab) ) then
      call display_string ( string, before='Declare ' )
      call dump_a_decl ( decl_table(num_decls), details=levels(tab) )
    end if
  end subroutine DECLARE

! -----------------------------------------------------  DECLARED  -----
  logical function DECLARED ( STRING )
    integer, intent(in) :: STRING
    if ( string > ubound(symbol_decl,1) ) & ! Assume string_table is increased
      & call increase_symbol_decl
    declared = symbol_decl(string) /= null_decl
  end function DECLARED

! ----------------------------------------------------  DUMP_DECL  -----
  subroutine DUMP_DECL ( Details )
    use Output_m, only: Output
    use String_Table, only: How_Many_Strings
    integer, intent(in), optional :: Details
    integer :: I    ! Loop inductor
    if ( num_decls <= 0 ) return
    call output ( ' dec  str', advance='yes' )
    do i = 1, how_many_strings()
      call dump_1_decl ( i, details=details )
    end do
  end subroutine DUMP_DECL

! --------------------------------------------------  DUMP_A_DECL  -----
  subroutine DUMP_A_DECL ( Decl, Before, Details )
    use Lexer_Core, only: Print_Source
    use Output_m, only: Output, Newline
    use String_Table, only: Display_String
    type(decls), intent(in) :: Decl
    character(len=*), intent(in), optional :: Before
    integer, intent(in), optional :: Details ! 0 -> no tree or source location
                                             ! >0 -> tree & source (default)
    integer :: MyDetails
    if ( present(before) ) call output ( before )
    myDetails = 1
    if ( present(details) ) myDetails = details
    if ( myDetails > 0 ) then
      call output ( ' type=' )
      call output ( trim(type_names(decl%type)) )
    end if
    select case ( decl%type )
    case ( action )
    case ( nonterminal )
    case ( terminal )
    case ( vocabulary )
    end select
    if ( decl%value /= 0 ) call display_string ( decl%value, before=' with value=' )
    if ( myDetails > 0 ) then
      if ( decl%tree /= null_tree ) then
        call output ( decl%tree, before=', tree=' )
        call print_source ( decl%tree, before=' ' )
      end if
    end if
    call newLine
  end subroutine DUMP_A_DECL

! --------------------------------------------------  DUMP_1_DECL  -----
  subroutine DUMP_1_DECL ( SYMBOL, Details )
    use Output_m, only: Output
    use String_Table, only: Display_String
    integer, intent(in) :: SYMBOL  ! Index of symbol whose declaration to dump
    integer, intent(in), optional :: Details
    integer :: DECL                ! Index of decl of "symbol"
    if ( symbol > ubound(symbol_decl,1) ) & ! Assume string_table is increased
      & call increase_symbol_decl
    decl = symbol_decl(symbol)
    if ( decl == null_decl ) then
      call output ( symbol, 5 ); call output ( ': ' )
      call display_string ( symbol )
      call output ( ' Not declared', advance='yes' )
    end if
    do while ( decl /= null_decl )
      call output ( decl, 4 )
      call output ( symbol, 5 ); call output ( ': ' )
      call display_string ( symbol )
      call dump_a_decl ( decl_table(decl), details=details )
      decl = decl_table(decl)%prior
    end do
  end subroutine DUMP_1_DECL

! -----------------------------------------------------  GET_DECL  -----
  type(decls) function Get_Decl ( String, Type ) result ( Decl )
  ! Get the latest declaration of "string" having a "type" field equal
  ! to "type".  Return it if any, else return the decl at null_decl.
    integer, intent(in) :: STRING  ! Index of string
    integer, intent(in) :: TYPE    ! "type" value to look for

    integer :: Where

    where = symbol_decl(string)
    do while ( where /= null_decl )
      decl = decl_table(where)
      if ( decl%type == type ) return
      where = decl_table(where)%prior
    end do
  end function Get_Decl

! ----------------------------------------------------  REDECLARE  -----
  subroutine REDECLARE ( STRING, TYPE, TREE, VALUE )
  ! Find the latest declaration for "string" of type "type".  If there
  ! isn't one, declare it.  Otherwise, change the "value", "units",
  ! "tree" and "values" fields of the found one.
    use String_Table, only: Display_String
    use Toggles, only: Tab, Toggle
    integer, intent(in) :: STRING  ! String index of name to declare
    integer, intent(in) :: TYPE    ! Type of object, e.g. UNITS, LABEL...
    integer, intent(in), optional :: TREE    ! Index of tree node of declaration
    integer, intent(in) :: VALUE   ! Declared value
    integer :: PRIOR
    if ( string > ubound(symbol_decl,1) ) & ! Assume string_table is increased
      & call increase_symbol_decl
    prior = symbol_decl(string)
    do
      if ( prior == null_decl ) then
        call declare ( string, type, tree, value )
        return
      end if
      if ( decl_table(prior)%type == type ) then
        decl_table(prior)%value = value
        if ( present(tree) )   decl_table(prior)%tree = tree
        if ( toggle(tab) ) then
          call display_string ( string, before='Redeclare ' )
          call dump_a_decl ( decl_table(prior), before=' with' )
        end if
        return
      end if
      prior = decl_table(prior)%prior
    end do
  end subroutine REDECLARE

! ----------------------------------------------------  INIT_DECL  -----
  subroutine INIT_DECL
                            !      type     value    tree   prior
    decl_table(null_decl) = decls( null_decl, 0, null_tree, null_decl )
    num_decls = 0
  end subroutine INIT_DECL

! =====     Private Procedures     =====================================
! -----------------------------------------  Increase_Symbol_Decl  -----
  subroutine Increase_Symbol_Decl ( Status )
  ! Double the size of the Symbol_Decl table.
    use, intrinsic :: ISO_Fortran_ENV, only: Error_Unit
    use String_Table, only: String_Table_Size
    integer, intent(out), optional :: Status
    integer :: N
    integer, allocatable :: New_Decl(:)
    integer :: Stat
    character(120) :: ERMSG ! From allocate

    n = ubound(symbol_decl,1)
    allocate ( new_decl(0:string_table_size()), stat=stat, errmsg=ermsg )
    if ( stat /= 0 ) then
      if ( present(status) ) then
        status = stat
        return
      end if
      write ( error_unit, '(a,i0/a)' ) &
        & 'DECL_TABLE%ALLOCATE_DECL-E- Unable to allocate NEW_DECL; STAT = ', &
        & stat, trim(ermsg)
      stop
    end if
    new_decl(0:n) = symbol_decl
    new_decl(n+1:) = null_decl
    call move_alloc ( new_decl, symbol_decl )
  end subroutine Increase_Symbol_Decl

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Declaration_Table

! $Log$
! Revision 1.1  2014/01/14 00:14:57  vsnyder
! Initial commit of new module for new LR
!
