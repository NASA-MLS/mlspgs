! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DECLARATION_TABLE

! It is necessary that the INIT_TABLES procedure from INIT_TABLES_MODULE
! is called before any procedures here.  The INIT_TABLES_MODULE depends on
! the application -- it is not part of the library -- but it is responsible
! for initializing the representation of intrinsic types accessed from
! the INTRINSIC module.

  use INTRINSIC, only: PHYQ_INDICES, PHYQ_INVALID, T_BOOLEAN, &
   &                   T_NUMERIC, T_NUMERIC_RANGE, T_STRING
  use MACHINE, only: IO_ERROR
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: DISPLAY_STRING, HOW_MANY_STRINGS, STRING_TABLE_SIZE
  use TOGGLES, only: TAB, TOGGLE

  implicit NONE
  private

  public :: ALLOCATE_DECL, DEALLOCATE_DECL, DECLARATION, DECLARE
  public :: DECLARED, DECLS, DUMP_DECL, DUMP_1_DECL, EMPTY
  public :: ENUM_VALUE, EXPRN, FIELD, GET_DECL, INIT_DECL, LABEL, LOG_VALUE
  public :: NAMED_VALUE, NULL_DECL, NUM_VALUE, PRIOR_DECL, RANGE, REDECLARE
  public :: SECTION, SECTION_NODE, STR_RANGE, STR_VALUE, SPEC, TREE_NODE
  public :: TYPE_MAP, TYPE_NAME, TYPE_NAMES, UNDECLARED, UNITS_NAME

  type :: DECLS
    double precision :: VALUE
    integer :: TYPE           ! "units", "variable", "spec", "field", "value", ...
    integer :: UNITS          ! Depends on "type" field:
                              ! "units_name" => Index of "units" of name,
                              !            e.g. km = length, ...
                              ! "spec" => Index of specification
                              ! "field" => Index of field
                              ! "section" => Index of section
    integer :: TREE           ! Index of declaration in the tree
    integer :: PRIOR          ! Index of previous declaration
  end type DECLS

  integer, parameter :: NULL_DECL = 0   ! Index and type of the null
                                        ! declaration

! Values of the "type" field of "decls":
  integer, parameter :: EMPTY = 0       ! The "type" field of the sentinel
  integer, parameter :: ENUM_VALUE = 1  ! An enumerator
  integer, parameter :: EXPRN = 2       ! The "tree" field points to an
                                        ! unevaluated expression
  integer, parameter :: FIELD = 3       ! Field of a structure definition
  integer, parameter :: LABEL = 4       ! A "name:" label for a stru
  integer, parameter :: LOG_VALUE = 5   ! Entity is a logical value, value is
                                        ! .false. if the "value" field is zero.
  integer, parameter :: NAMED_VALUE = 6 ! X = expr
  integer, parameter :: NUM_VALUE = 7   ! Entity is a numeric value, value is
                                        ! in the "value" field"
  integer, parameter :: RANGE = 8       ! A range -- not used in decl table
  integer, parameter :: SECTION = 9     ! Name of a section
  integer, parameter :: SECTION_NODE = 10 ! Tree node of a section
  integer, parameter :: STR_RANGE = 11  ! String range -- for dates
  integer, parameter :: STR_VALUE = 12  ! The string is the value
  integer, parameter :: SPEC = 13       ! Name of a specification, e.g. vGrid
  integer, parameter :: TREE_NODE = 14  ! Name of a tree node, e.g. n_plus
  integer, parameter :: TYPE_NAME = 15  ! Name of a data type
  integer, parameter :: UNDECLARED = 16 ! Entity is undeclared
  integer, parameter :: UNITS_NAME = 17 ! Name is a units name, e.g. km, hPa
                                        ! Scale to "canonical" units of the
                                        ! name is in "value", e.g. km = 1000.0

  character(len=*), parameter :: TYPE_NAMES(empty:units_name) = &
  (/ 'empty     ', 'enum_value', 'exprn     ', 'field     ', 'label     ', &
     'log_value ', 'nam_value ', 'num_value ', 'range     ', 'section   ', &
     'section_n ', 'str_range ', 'str_value ', 'spec      ', 'tree      ', &
     'type_name ', 'undeclared', 'units     ' /)

! Mapping from declaration table types to data types:
  integer, parameter :: TYPE_MAP(empty:units_name) = &
  (/ 0, 0, 0, 0, 0, t_boolean, 0, t_numeric, t_numeric_range, 0, 0, 0, &
     t_string, 0, 0, 0, 0, 0 /)

! -----     Private declarations     -----------------------------------
  type(decls), save, allocatable :: DECL_TABLE(:)
  integer, save :: NUM_DECLS            ! amount of decl_table used
  integer, save, allocatable :: SYMBOL_DECL(:) ! indexed by string
                                        ! index, gives index in decl_table.

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================
! ------------------------------------------------  ALLOCATE_DECL  -----
  subroutine ALLOCATE_DECL ( NDECLS, STAT )
  ! Allocate NDECLS declarations.  Allocate STRING_TABLE_SIZE symbol
  ! declarations (indexes from STRING_TABLE to DECL_TABLE).  Also does
  ! INIT_DECL.
    integer, intent(in) :: NDECLS       ! Number to allocate
    integer, intent(out), optional :: STAT   ! From ALLOCATE statement
    integer :: MY_STAT
    if ( allocated(decl_table) ) then; deallocate(decl_table); end if
    if ( allocated(symbol_decl) ) then; deallocate(symbol_decl); end if
    allocate ( decl_table(0:ndecls), stat=my_stat )
    if ( my_stat /= 0 ) then
      if ( present(stat) ) then
        stat = my_stat
        return
      end if
      call io_error &
      ( 'DECL_TABLE%ALLOCATE_DECL-E- Unable to allocate storage', my_stat )
      stop
    end if
    allocate ( symbol_decl(0:string_table_size()), stat=my_stat )
    if ( my_stat /= 0 ) then
      if ( present(stat) ) then
        stat = my_stat
        return
      end if
      call io_error &
      ( 'DECL_TABLE%ALLOCATE_DECL-E- Unable to allocate storage', my_stat )
      stop
    end if
    call init_decl
    symbol_decl = null_decl
    return
  end subroutine ALLOCATE_DECL
! ----------------------------------------------  DEALLOCATE_DECL  -----
  subroutine DEALLOCATE_DECL
    if ( allocated(decl_table) ) deallocate ( decl_table )
    if ( allocated(symbol_decl) ) deallocate( symbol_decl )
  end subroutine DEALLOCATE_DECL
! --------------------------------------------------  DECLARATION  -----
  type(decls) function DECLARATION ( STRING )
    integer, intent(in) :: STRING  ! String index for which declaration needed
    if ( string > ubound(symbol_decl,1) ) & ! Assume string_table is increased
      & call increase_symbol_decl
    declaration = decl_table(symbol_decl(string))
  end function DECLARATION
! ------------------------------------------------------  DECLARE  -----
  subroutine DECLARE ( STRING, VALUE, TYPE, UNITS, TREE )
    integer, intent(in) :: STRING  ! String index of name to declare
    double precision, intent(in) :: VALUE    ! Declared value
    integer, intent(in) :: TYPE    ! Type of object, e.g. UNITS, LABEL...
    integer, intent(in) :: UNITS   ! Units of value -- index of units string
    integer, intent(in) :: TREE    ! Index of tree node of declaration

    integer :: STAT
    type(decls), allocatable :: OLD_DECL(:)

    num_decls = num_decls + 1
    if ( num_decls > ubound(decl_table,1) ) then
    ! Double size of declaration table
      allocate ( old_decl(0:ubound(decl_table,1)), stat=stat )
      if ( stat /= 0 ) then
        call io_error &
        ( 'DECLARATION_TABLE%DECLARE-E- Unable to allocate storage', stat )
        stop
      end if
      old_decl = decl_table
      deallocate ( decl_table )
      allocate( decl_table(0:2*size(old_decl)), stat=stat )
      if ( stat /= 0 ) then
        call io_error &
        ( 'DECLARATION_TABLE%DECLARE-E- Unable to allocate storage', stat )
        stop
      end if
      decl_table(0:ubound(old_decl,1)) = old_decl
      deallocate ( old_decl )
    end if
    if ( string > ubound(symbol_decl,1) ) & ! Assume string_table is increased
      & call increase_symbol_decl
    decl_table(num_decls) = decls ( value, type, units, tree, &
                                    symbol_decl(string) )
    symbol_decl(string) = num_decls
    if ( toggle(tab) ) then
      call output ( 'Declare ' ); call display_string ( string )
      call output ( ' with VALUE = ' ); call output ( value )
      call output ( ', TYPE = ' ); call output ( type )
      call output ( ', UNITS = ' ); call output ( units )
      call output ( ', TREE = ' ); call output ( tree, advance='yes' )
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
  subroutine DUMP_DECL
    integer :: I    ! Loop inductor
    if ( num_decls <= 0 ) return
    call output ( ' dec str', advance='yes' )
    do i = 1, how_many_strings()
      call dump_1_decl ( i )
    end do
  end subroutine DUMP_DECL
! --------------------------------------------------  DUMP_1_DECL  -----
  subroutine DUMP_1_DECL ( SYMBOL )
    integer, intent(in) :: SYMBOL  ! Index of symbol whose declaration to dump
    integer :: DECL                ! Index of decl of "symbol"
    if ( symbol > ubound(symbol_decl,1) ) & ! Assume string_table is increased
      & call increase_symbol_decl
    decl = symbol_decl(symbol)
    do while ( decl /= null_decl )
      call output ( decl, 4 )
      call output ( symbol, 4 ); call output ( ': ' )
      call display_string ( symbol )
      call output ( ' value=' )
      call output ( decl_table(decl)%value )
      call output ( ' type=' )
      call output ( trim(type_names(decl_table(decl)%type)) )
      call output ( ' units=' )
      if ( decl_table(decl)%type == units_name ) then
        call display_string ( phyq_indices(decl_table(decl)%units) )
      else
        call output ( decl_table(decl)%units )
      end if
      call output ( ' tree=' )
      call output ( decl_table(decl)%tree, advance='yes' )
      decl = decl_table(decl)%prior
    end do
  end subroutine DUMP_1_DECL
! -----------------------------------------------------  GET_DECL  -----
  type(decls) function GET_DECL ( STRING, TYPE, UNITS, TREE )
  ! Get the latest declaration of "string" having a "type" field equal
  ! to "type" (if present), a "units" field equal to "units" (if
  ! present) and a "tree" field equal to "tree" (if present).  Return it
  ! if any, else return the decl at null_decl.
    integer, intent(in) :: STRING  ! Index of string
    integer, intent(in), optional :: TYPE  ! "type" value to look for
    integer, intent(in), optional :: UNITS ! "units" value to look for
    integer, intent(in), optional :: TREE  ! "tree" value to look for

    logical :: GOT_IT
    integer :: PRIOR

    get_decl = declaration(string)
    do
      got_it = .true.
      if ( present(type) ) got_it = get_decl%type == type
      if ( got_it .and. present(units) ) got_it = get_decl%units == units
      if ( got_it .and. present(tree) ) got_it = get_decl%tree == tree
      if ( got_it ) return
      prior = get_decl%prior
      get_decl = prior_decl(get_decl)
      if ( prior == null_decl ) return
    end do
  end function GET_DECL
! ----------------------------------------------------  INIT_DECL  -----
  subroutine INIT_DECL
    !                             value type  units     tree prior
    decl_table(null_decl) = decls(0.0d0,empty,phyq_invalid,0,null_decl)
    num_decls = 0
  end subroutine INIT_DECL
! ---------------------------------------------------  PRIOR_DECL  -----
  type(decls) function PRIOR_DECL ( THE_DECL, TYPE, UNITS )
  ! Return the prior declaration of "the_decl" with specified values of
  ! "type" and "units" (if they're specified).
    type(decls), intent(in) :: THE_DECL
    integer, intent(in), optional :: TYPE, UNITS
    integer :: PRIOR
    prior = the_decl%prior
    do
      prior_decl = decl_table(prior)
      if ( prior_decl%type == empty ) return
      if ( present(type) ) then
        if ( prior_decl%type == type ) then
          if ( .not. present(units) ) return
          if ( prior_decl%units == units ) return
        end if
      else
        if ( .not. present(units) ) return
        if ( prior_decl%units == units ) return
      end if
      prior = prior_decl%prior
    end do
  end function PRIOR_DECL
! ----------------------------------------------------  REDECLARE  -----
  subroutine REDECLARE ( STRING, VALUE, TYPE, UNITS, TREE )
  ! Find the latest declaration for "string" of type "type".  If there
  ! isn't one, declare it.  Otherwise, change the "value", "units" and
  ! "tree" fields of the found one.
    integer, intent(in) :: STRING  ! String index of name to declare
    double precision, intent(in) :: VALUE    ! Declared value
    integer, intent(in) :: TYPE    ! Type of object, e.g. UNITS, LABEL...
    integer, intent(in) :: UNITS   ! Units of value -- index of units string
    integer, intent(in) :: TREE    ! Index of tree node of declaration
    integer :: PRIOR
    if ( string > ubound(symbol_decl,1) ) & ! Assume string_table is increased
      & call increase_symbol_decl
    prior = symbol_decl(string)
    do
      if ( prior == null_decl ) then
        call declare ( string, value, type, units, tree )
        return
      end if
      if ( decl_table(prior)%type == type ) then
        decl_table(prior)%value = value
        decl_table(prior)%units = units
        decl_table(prior)%tree = tree
        return
      end if
      prior = decl_table(prior)%prior
    end do
  end subroutine REDECLARE

! =====     Private Procedures     =====================================
! -----------------------------------------  Increase_Symbol_Decl  -----
  subroutine Increase_Symbol_Decl ( Status )
  ! Double the size of the Symbol_Decl table.
    integer, intent(out), optional :: Status
    integer, allocatable :: Old_Decl(:)
    integer :: Stat

    allocate ( old_decl(0:ubound(symbol_decl,1)), stat=stat )
    if ( stat /= 0 ) then
      if ( present(status) ) then
        status = stat
        return
      end if
      call io_error &
      ( 'DECL_TABLE%Increase_Symbol_Decl-E- Unable to allocate storage', stat )
      stop
    end if
    old_decl = symbol_decl
    deallocate ( symbol_decl )
    allocate ( symbol_decl(0:string_table_size()), stat=stat )
    if ( stat /= 0 ) then
      if ( present(status) ) then
        status = stat
        return
      end if
      call io_error &
      ( 'DECL_TABLE%Increase_Symbol_Decl-E- Unable to allocate storage', stat )
      stop
    end if
    symbol_decl(0:ubound(old_decl,1)) = old_decl
    symbol_decl(ubound(old_decl,1)+1:) = null_decl
    deallocate ( old_decl )
  end subroutine Increase_Symbol_Decl
end module DECLARATION_TABLE

! $Log$
! Revision 2.3  2001/04/16 23:05:29  vsnyder
! SAVE some module variables
!
! Revision 2.2  2001/04/05 01:28:06  vsnyder
! Add 'increase the symbol declaration table size' code
!
! Revision 2.1  2000/10/11 18:57:28  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:49  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
