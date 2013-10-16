! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module DECLARATION_TABLE

! It is necessary that the INIT_TABLES procedure from INIT_TABLES_MODULE
! is called before any procedures here.  The INIT_TABLES_MODULE depends on
! the application -- it is not part of the library -- but it is responsible
! for initializing the representation of intrinsic types accessed from
! the INTRINSIC module.

  use INTRINSIC, only: LIT_INDICES, PHYQ_INDICES, PHYQ_INVALID, T_A_DOT_B, &
    &                  T_BOOLEAN, T_NUMERIC, T_NUMERIC_RANGE, T_STRING
  use MACHINE, only: IO_ERROR
  use OUTPUT_M, only: NewLine, OUTPUT
  use STRING_TABLE, only: CREATE_STRING, DISPLAY_STRING, HOW_MANY_STRINGS, &
    &                     STRING_TABLE_SIZE
  use TOGGLES, only: TAB, TOGGLE
  use Tree, only: Null_Tree

  implicit NONE
  private

  public :: ALLOCATE_DECL, Allocate_Test, Value_Allocate, DEALLOCATE_DECL
  public :: Deallocate_Test, DECLARATION, DECLARE, DECLARED, DECLS, DOT
  public :: DUMP_DECL, DUMP_A_DECL, DUMP_1_DECL, DUMP_VALUES, EMPTY, ENUM_VALUE
  public :: EXPRN, EXPRN_M, EXPRN_V, FIELD, FUNCTION, GET_DECL, GET_TYPE
  public :: INIT_DECL, LABEL, LOG_VALUE, NAMED_VALUE, NULL_DECL, NUM_VALUE
  public :: PHYS_UNIT_NAME, PRIOR_DECL, RANGE, REDECLARE, SECTION
  public :: SECTION_NODE, STR_RANGE, STR_VALUE, SPEC, TREE_NODE, TYPE_MAP
  public :: TYPE_NAME, TYPE_NAMES, TYPE_NAME_INDICES, UNDECLARED, UNITS_NAME
  public :: Value_Add, Value_Deallocate, Value_Diff, Value_Div, Value_Equal
  public :: Value_Neg, Value_Prod, Value_Same, Value_t, Value_Unequal, VARIABLE

  public :: Operator(==), Operator(/=)
  public :: Operator(+), Operator(-), Operator(*), Operator(/)

  include "Value_T_Interfaces.f9h"

  interface Allocate_Test
    module procedure Value_Allocate
  end interface

  interface Deallocate_Test
    module procedure Value_Deallocate
  end interface

  integer, parameter :: NULL_DECL = 0   ! Index and type of the null
                                        ! declaration

! Values of the "type" field of "decls":
  integer, parameter :: EMPTY = 0       ! The "type" field of the sentinel
  integer, parameter :: DOT = 1         ! A.B -- not used in decl table
  integer, parameter :: ENUM_VALUE = 2  ! An enumerator
  integer, parameter :: EXPRN = 3       ! The "tree" field points to a
                                        ! scalar expression
  integer, parameter :: EXPRN_M = 4     ! The "tree" field points to a
                                        ! matrix expression
  integer, parameter :: EXPRN_V = 5     ! The "tree" field points to a
                                        ! vector expression
  integer, parameter :: FIELD = 6       ! Field of a structure definition
  integer, parameter :: FUNCTION = 7    ! Name is a built-in function
  integer, parameter :: LABEL = 8       ! A "name:" label for a stru
  integer, parameter :: LOG_VALUE = 9   ! Entity is a logical value, value is
                                        ! .false. if the "value" field is zero.
  integer, parameter :: NAMED_VALUE = 10! X = expr
  integer, parameter :: NUM_VALUE = 11  ! Entity is a numeric value, value is
                                        ! in the "value" field"
  integer, parameter :: PHYS_UNIT_NAME = 12 ! PHYQ_....
  integer, parameter :: RANGE = 13      ! A range -- not used in decl table
  integer, parameter :: SECTION = 14    ! Name of a section
  integer, parameter :: SECTION_NODE = 15 ! Tree node of a section
  integer, parameter :: STR_RANGE = 16  ! String range -- for dates
  integer, parameter :: STR_VALUE = 17  ! The string is the value
  integer, parameter :: SPEC = 18       ! Name of a specification, e.g. vGrid
  integer, parameter :: TREE_NODE = 19  ! Name of a tree node, e.g. n_plus
  integer, parameter :: TYPE_NAME = 20  ! Name of a data type
  integer, parameter :: UNDECLARED = 21 ! Entity is undeclared
  integer, parameter :: UNITS_NAME = 22 ! Name is a units name, e.g. km, hPa
                                        ! Scale to "canonical" units of the
                                        ! name is in "value", e.g. km = 1000.0
  integer, parameter :: VARIABLE = 23   ! Name is a variable, e.g. A := <expr>

  integer, parameter :: LAST_TYPE = VARIABLE

  type :: Value_t             ! Used if decls%type == variable or named_value
    integer :: TYPE = empty   ! Value type, subset of decls%type
    double precision :: VALUE(2) = 0.0d0 ! Real if type == num_value or range
                              ! 0 => false, 1 => true if type == log_value
                              ! string index if type == enum_value, label,
                              !   named_value, str_range, str_value
    integer :: Units(2) = 0   ! Subset of decls%units, depends on Type
  end type Value_t

  ! Type            Units                Value           Tree
  ! EMPTY           0                    0.0d0           null_tree
  ! DOT
  ! ENUM_VALUE      data type index      string index    DT_Def
  ! EXPRN
  ! EXPRN_M
  ! EXPRN_V
  ! FIELD           field index          0.0d0           Field_type
  !                                      0.0d0           Field_Spec
  ! FUNCTION        function index       string index    Func_def
  ! LABEL           PHYQ_Invalid         0.0d0           Spec_Args
  ! LOG_VALUE
  ! NAMED_VALUE     param index          string index    Name_Def
  ! NUM_VALUE       PHYQ_Dimensionless   its value       Number
  ! PHYS_UNIT_NAME  its lit index        value           string index
  ! RANGE
  ! SECTION         section index        string index    Section
  ! SECTION_NODE
  ! STR_RANGE
  ! STR_VALUE       PHYQ_Invalid         0.0d0           String
  ! SPEC            spec index           string index    Spec_def
  ! TREE_NODE
  ! TYPE_NAME       type index           string index    DT_Def
  ! UNDECLARED      PHYQ_Invalid         0.0d0           Identifier
  ! UNITS_NAME      PHYQ_...             Value           string index
  ! VARIABLE        Type of first        string index    Identifier
  !                 element of value

  type :: DECLS
    double precision :: VALUE = 0.0d0
    integer :: TYPE = empty   ! "units", "variable", "spec", "field", "value", ...
    integer :: UNITS = 0      ! Depends on "type" field:
                              ! "units_name" => Index of "units" of name,
                              !            e.g. km = length, ...
                              ! "enum_value" => The enumerator's lit index
                              ! "exprn" => < 0 offset in matrix database
                              !            > 0 offset in vector database
                              !            0 = double-precision value only
                              ! "field" => Index of field
                              ! "function" => Index of function
                              ! "named_value" => Type of its first value
                              ! "section" => Index of section
                              ! "spec" => Index of specification
                              ! "variable" => Type of its first value
    integer :: TREE = null_tree ! Index of declaration in the tree, except for
                              ! "units_name" => string index of name
    integer :: PRIOR = null_decl ! Index of previous declaration
    type(value_t), allocatable :: Values(:) ! If type == variable or named_value
  end type DECLS

  character(len=*), parameter :: TYPE_NAMES(empty:variable) = &
  (/ 'empty         ', 'dot           ', 'enum_value    ', 'exprn         ', &
     'exprn_m       ', 'exprn_v       ', 'field         ', 'function      ', &
     'label         ', 'log_value     ', 'nam_value     ', 'num_value     ', &
     'phys_unit_name', 'range         ', 'section       ', 'section_n     ', &
     'str_range     ', 'str_value     ', 'spec          ', 'tree          ', &
     'type_name     ', 'undeclared    ', 'units         ', 'variable      ' /)

  ! String indices of type names, filled on first call to Get_Type or
  ! Init_Decl
  integer :: TYPE_NAME_INDICES(empty:variable) = -1

! Mapping from declaration table types to data types:
  integer, parameter :: TYPE_MAP(empty:variable) =        &
  (/ 0, t_a_dot_b, 0,         0,               0,         &
     0, 0,         0,         0,               t_boolean, &
     0, t_numeric, 0,         t_numeric_range, 0,         &
     0, 0,         t_string,  0,               0,         &
     0, 0,         0,         0 /)

! -----     Private declarations     -----------------------------------
  type(decls), save, allocatable :: DECL_TABLE(:)
  integer, save :: NUM_DECLS            ! amount of decl_table used
  integer, save, allocatable :: SYMBOL_DECL(:) ! indexed by string
                                        ! index, gives index in decl_table.

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
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
  subroutine DECLARE ( STRING, VALUE, TYPE, UNITS, TREE, VALUES )
    integer, intent(in) :: STRING  ! String index of name to declare
    double precision, intent(in) :: VALUE    ! Declared value
    integer, intent(in) :: TYPE    ! Type of object, e.g. UNITS, LABEL...
    integer, intent(in) :: UNITS   ! Units of value -- index of units string
    integer, intent(in) :: TREE    ! Index of tree node of declaration
    type(value_t), allocatable, intent(inout), optional :: Values(:)

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
                                    symbol_decl(string), null() )
    if ( present(values) ) call move_alloc ( values, decl_table(num_decls)%values )
    symbol_decl(string) = num_decls
    if ( toggle(tab) ) then
      call display_string ( string, before='Declare ' )
      call dump_a_decl ( decl_table(num_decls), before=' with' )
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
    call output ( ' dec  str', advance='yes' )
    do i = 1, how_many_strings()
      call dump_1_decl ( i )
    end do
  end subroutine DUMP_DECL

! --------------------------------------------------  DUMP_A_DECL  -----
  subroutine DUMP_A_DECL ( Decl, Before, Value_Only )
    use Lexer_Core, only: Print_Source
    use Tree, only: Subtree, Sub_Rosa
    type(decls), intent(in) :: Decl
    character(len=*), intent(in), optional :: Before
    logical, intent(in), optional :: Value_Only
    logical :: All
    all = .true.
    if ( present(value_only) ) all = .not. value_only
    if ( present(before) ) call output ( before )
    if ( all ) then
      call output ( ' type=' )
      call output ( trim(type_names(decl%type)) )
    end if
    if ( .not. allocated(decl%values) ) &
      & call output ( decl%value, before=' value=' )
    select case ( decl%type )
    case ( enum_value )
      call display_string ( sub_rosa(subtree(1,decl%tree)), before=' of type ' )
    case ( phys_unit_name )
      call display_string ( lit_indices(decl%units), before=' units=' )
    case ( named_value, variable )
      if ( decl%units >= empty .and. decl%units <= last_type ) &
        & call output ( ' ' // trim(type_names(decl%units)) )
      if ( decl%units == str_value .and. decl%value > 0 ) &
        & call display_string ( nint(decl%value), before=' ' )
    case ( units_name )
      if ( phyq_indices(decl%units) == 0 ) then
        call output ( ' units=<unknown> ' )
      else
        call display_string ( phyq_indices(decl%units), before=' units=' )
      end if
    case default
      call output ( decl%units, before=' units=' )
    end select
    if ( all ) then
      call output ( ' tree=' )
      call output ( decl%tree )
      call print_source ( decl%tree, before=' ' )
    end if
    if ( allocated(decl%values) ) then
      call output ( ' values:', advance='yes' )
      call dump_values ( decl%values )
    else
      call newLine
    end if
  end subroutine DUMP_A_DECL

! --------------------------------------------------  DUMP_1_DECL  -----
  subroutine DUMP_1_DECL ( SYMBOL )
    integer, intent(in) :: SYMBOL  ! Index of symbol whose declaration to dump
    integer :: DECL                ! Index of decl of "symbol"
    if ( symbol > ubound(symbol_decl,1) ) & ! Assume string_table is increased
      & call increase_symbol_decl
    decl = symbol_decl(symbol)
    do while ( decl /= null_decl )
      call output ( decl, 4 )
      call output ( symbol, 5 ); call output ( ': ' )
      call display_string ( symbol )
      call dump_a_decl ( decl_table(decl) )
      decl = decl_table(decl)%prior
    end do
  end subroutine DUMP_1_DECL

! --------------------------------------------------  DUMP_VALUES  -----
  subroutine DUMP_VALUES ( VALUES, BEFORE, ADVANCE )
    type(value_t), intent(in) :: Values(:)
    character(len=*), intent(in), optional :: BEFORE, ADVANCE
    integer :: I, N
    if ( present(before) ) call output ( before, advance )
    do i = 1, size(values)
      if ( size(values) > 1 ) then
        call output ( i, places=3 )
        call output ( ': ' )
      end if
      select case ( values(i)%type )
      case ( enum_value )  ! values%value(1) is a string index
        call display_string ( nint(values(i)%value(1)) , advance='yes' )
      case ( label )       ! values%value(1) is a string index
        call display_string ( nint(values(i)%value(1)), advance='yes' )
      case ( log_value )   ! values%value(1) 0= 0 => false, else true
        call output ( trim(merge('true ','false',values(i)%value(1)/=0.0d0)), &
          & advance='yes' )
      case ( named_value, variable )
        select case ( values(i)%units(1) )
        case ( num_value )
          call output ( values(i)%value(1), advance='yes' )
        case ( range )
          call output ( values(i)%value(1) )
          call output ( values(i)%value(2), before=' : ', advance='yes' )
        case ( str_range )
          n = nint(values(i)%value(1))
          if ( n > 0 ) call display_string ( nint(values(i)%value(1)) )
          n = nint(values(i)%value(2))
          if ( n > 0 ) then
            call display_string ( n, before=' : ', advance='yes' )
          else
            call output ( ' : ', advance='yes' )
          end if
        case ( str_value )
          n = nint(values(i)%value(1))
          if ( n > 0 ) then
            call display_string ( n, advance='yes' )
          else
            call newLine
          end if
        case default
          call newLine
        end select
      case ( num_value )   ! values%value(1) is a number
        call output ( values(i)%value(1), advance='yes' )
      case ( range )       ! values%value(1:2) are numbers
        call output ( values(i)%value(1) )
        call output ( values(i)%value(2), before=' : ', advance='yes' )
      case ( str_range )   ! values%value(1:2) are string indices
        call display_string ( nint(values(i)%value(1)) )
        call display_string ( nint(values(i)%value(2)), before=' : ', &
          & advance='yes' )
      case ( str_value )   ! values%value(1) is a string index
        call display_string ( nint(values(i)%value(1)), advance='yes' )
      case default
        call output ( trim(type_names(values(i)%type )), advance='yes' )
      end select
    end do
  end subroutine DUMP_VALUES

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

  ! --------------------------------------------------  Get_Type  -----
  ! Return the string index for the type indexed by Decor
  integer function Get_Type ( Decor )
    integer, intent(in) :: Decor  ! Tree node decoration
    if ( (empty) < 0 ) call init_type_indices
    get_type = 0
    if ( decor >= empty .and. decor <= last_type ) &
      & get_type = type_name_indices ( decor )
  end function Get_Type

! ----------------------------------------------------  INIT_DECL  -----
  subroutine INIT_DECL
    !                             value type  units     tree prior
    decl_table(null_decl) = decls(0.0d0,empty,phyq_invalid,0,null_decl)
    num_decls = 0
    if ( type_name_indices(empty) < 0 ) call init_type_indices
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
  subroutine REDECLARE ( STRING, VALUE, TYPE, UNITS, TREE, VALUES )
  ! Find the latest declaration for "string" of type "type".  If there
  ! isn't one, declare it.  Otherwise, change the "value", "units",
  ! "tree" and "values" fields of the found one.
    integer, intent(in) :: STRING  ! String index of name to declare
    double precision, intent(in) :: VALUE    ! Declared value
    integer, intent(in) :: TYPE    ! Type of object, e.g. UNITS, LABEL...
    integer, intent(in), optional :: UNITS   ! Units of value -- index of units string
    integer, intent(in), optional :: TREE    ! Index of tree node of declaration
    type(value_t), allocatable, intent(inout), optional :: Values(:)
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
        if ( present(units) )  decl_table(prior)%units = units
        if ( present(tree) )   decl_table(prior)%tree = tree
        if ( present(values) ) then
          call deallocate_test ( decl_table(prior)%values, &
            & 'decl_table(prior)%values', moduleName )
          call move_alloc ( values, decl_table(prior)%values )
        end if
        if ( toggle(tab) ) then
          call display_string ( string, before='Redeclare ' )
          call dump_a_decl ( decl_table(prior), before=' with' )
        end if
        return
      end if
      prior = decl_table(prior)%prior
    end do
  end subroutine REDECLARE

! -----------------------------------------------  Value_Allocate  -----
  subroutine Value_Allocate ( Value, N, ItsName, ModuleName )
    use Allocate_Deallocate, only: Test_Allocate
    type(value_t), allocatable :: Value(:)
    integer, intent(in) :: N
    character(len=*), intent(in) :: ItsName, ModuleName
    integer :: Stat
    call deallocate_test ( value, itsName, moduleName )
    allocate ( value(1:n), stat=stat )
    call test_allocate ( stat, moduleName, ItsName, [1], [n], &
      & storage_size(value) / 8 )
  end subroutine Value_Allocate

! ---------------------------------------------  Value_Deallocate  -----
  subroutine Value_Deallocate ( Value, ItsName, ModuleName )
    use Allocate_Deallocate, only: Memory_Units, Test_Deallocate
    type(value_t), allocatable :: Value(:)
    character(len=*), intent(in) :: ItsName, ModuleName
    integer :: N, Stat
    if ( allocated(value) ) then
      n = size(value)
      deallocate ( value, stat=stat )
      call test_deallocate ( stat, moduleName, ItsName,  &
        & n * (storage_size(value) / 8) / memory_units )
    end if
  end subroutine Value_Deallocate

  include "Value_T_Implementations.f9h"

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
! --------------------------------------------  INIT_TYPE_INDICES  -----
  subroutine INIT_TYPE_INDICES
    integer :: I
    do i = empty, last_type
      type_name_indices(i) = create_string ( trim(type_names(i)) )
    end do
  end subroutine INIT_TYPE_INDICES

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module DECLARATION_TABLE

! $Log$
! Revision 2.16  2013/10/16 01:12:47  vsnyder
! Cannonball polishing
!
! Revision 2.15  2013/10/11 00:44:28  vsnyder
! Get Values_t operations from the include files
!
! Revision 2.14  2013/10/09 23:39:49  vsnyder
! Add Allocate_Test etc, equality etc.
!
! Revision 2.13  2013/10/02 01:30:03  vsnyder
! Add 'variable' type
!
! Revision 2.12  2013/09/19 23:29:33  vsnyder
! Add Dump_A_Decl, PHYS_Unit
!
! Revision 2.11  2013/09/17 00:55:07  vsnyder
! Add Dot type
!
! Revision 2.10  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.9  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.8  2004/01/23 05:35:50  livesey
! Made dump_1_decl a little more robust
!
! Revision 2.7  2004/01/17 03:04:48  vsnyder
! Provide for functions in expressions
!
! Revision 2.6  2004/01/16 23:51:23  vsnyder
! Add more declaration table types for Algebra
!
! Revision 2.5  2002/10/08 00:09:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2002/07/18 22:04:06  vsnyder
! Improve some debugging print
!
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
