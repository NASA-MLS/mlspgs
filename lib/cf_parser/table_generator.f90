module TABLE_GENERATOR

! Traverse the tree output by the parser, as annotated by tree_checker.
! Output tables in the "mlscf" format.

  use DECLARATION_TABLE, only: DUMP_1_DECL, EMPTY, ENUM_VALUE, LABEL, &
    & LOG_VALUE, RANGE, STR_RANGE, STR_VALUE, UNDECLARED
  use EXPR_M, only: EXPR
  use INIT_TABLES_MODULE ! everything beginning with PHYQ_
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCF, only: ALLOCATE_MLSCF, MLSCF_T, MLSCFCELL_T
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: GET_STRING
  use TRACE_M, only: DEPTH, TRACE_BEGIN, TRACE_END
  use TOGGLES, only: GEN, TOGGLE
  use TREE, only: DUMP_TREE_NODE, NODE_ID, &
                  NODE_KIND, NSONS, PSEUDO, SOURCE_REF, &
                  SUB_ROSA, SUBTREE
  use TREE_TYPES ! Everything, especially everything beginning with N_
  implicit NONE
  private

  public :: GENERATE_TABLE

! -----     Private declarations     -----------------------------------
  integer, private :: ERROR   ! 0 => No errors

  ! Counters for filling the tables
  integer, private :: ISECTION     ! Section number
  integer, private :: IDEF         ! Definition (X=...) in the section
  integer, private :: ISPEC        ! Specification (Y,...) number in the section
  integer, private :: ICELL        ! Cell (field=expr) number in the spec

! Error codes for "announce_error"
  integer, private, parameter :: NO_CODE_FOR = 1
  integer, private, parameter :: TOO_MANY_CELLS = 2
  integer, private, parameter :: TOO_MANY_DEFS = 3
  integer, private, parameter :: TOO_MANY_SPECS = 4

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================
! -----------------------------------------------  GENERATE_TABLE  -----
  subroutine GENERATE_TABLE ( ROOT, HOW_MANY_SECTIONS, L2CF_DATA, UPCASE )
  ! Traverse the abstract syntax tree starting at ROOT, which should
  ! be a N_CFS node (but we don't check).
    integer, intent(in) :: ROOT               ! Root of tree output by parser
                                              ! and decorated by check_tree
    integer, intent(in) :: HOW_MANY_SECTIONS  ! How many begin-end blocks
                                              ! were found by check_tree
    type(mlscf_t), intent(inout) :: L2CF_DATA ! Tables
    logical, intent(in), optional :: UPCASE   ! Upper-case strings if present
                                              ! and true.

    integer :: I              ! Loop inductor
    logical :: MYUPCASE
    integer :: SON            ! Son of root

    depth = 0
    error = 0
    isection = 0
    if ( toggle(gen) ) call trace_begin ( 'GENERATE_TABLE', root )
    myupcase = .false.
    if ( present(upcase) ) myupcase = upcase
    nullify ( l2cf_data%sections )
    call allocate_mlscf ( l2cf_data, how_many_sections )
    do i = 1, nsons(root)
      son = subtree(i,root)
      select case ( node_id(son) )
      case ( n_cf );        call one_cf ( son )
      case ( n_dt_def );    ! nothing needed
      case ( n_section );   ! nothing needed
      case ( n_spec_def );  ! nothing needed
      case default ;        call announce_error ( son, no_code_for )
      end select
    end do
    if ( toggle(gen) ) call trace_end ( 'GENERATE_TABLE' )
! =====     Internal Procedures     ====================================
  contains
! -----------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE, SONS )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message
    integer, intent(in), optional :: SONS(:) ! Tree nodes, maybe sons of
                                   ! "where".  If they're pseudo_terminal,
                                   ! their declarations are dumped.
    integer :: I                   ! Index for "sons"

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref(where) )
    call output ( ': ' )
    select case ( code )
    case ( no_code_for )
      call output ( ' there is no code to analyze ' )
      call dump_tree_node ( where, 0, advance='yes' )
    case ( too_many_cells )
      call output ( '.  Too many cells.', advance = 'yes' )
    case ( too_many_defs )
      call output ( '.  Too many defs.', advance = 'yes' )
    case ( too_many_specs )
      call output ( '.  Too many specs.', advance = 'yes' )
    case default
      call output ( 'No message in TABLE_GENERATOR for error code ' )
      call output ( code, advance='yes' )
      stop
    end select
    if ( present(sons) ) then
      do i = 1, size(sons)
        if ( node_kind(sons(i)) == pseudo ) &
          call dump_1_decl ( sub_rosa(sons(i)) )
      end do
    end if
  end subroutine ANNOUNCE_ERROR
! ----------------------------------------------------  FILL_CELL  -----
  subroutine FILL_CELL ( ROOT, CELL, IC, MSGNUM )
  ! Fill a cell from "name = expr", "/name" or "expr"
    integer, intent(in) :: ROOT              ! Root of subtree
    type(MlscfCell_T), intent(out) :: CELL(:)  ! Cell array
    integer, intent(inout) :: IC             ! Index of CELL, updated here
    integer, intent(in) :: MSGNUM            ! Message if IC > MAXIC

    integer :: GSON1, GSON2   ! Sons of Son
    integer :: I              ! Index of son
    integer :: SON            ! A son of root
    integer :: TYPE           ! Type of son
    integer :: UNITS(2)       ! Units of son
    double precision :: VALUE(2) ! Value of son

    if ( toggle(gen) ) call trace_begin ( 'FILL_CELL', root )
    ic = ic + 1
    if ( ic > size(cell) ) then
      call announce_error ( root, msgnum )
      call output ( 'Number of them = ' ); call output ( ic, advance='yes' )
      stop
    end if
    cell(ic)%charvalue = ' '
    cell(ic)%units = 0
    cell(ic)%rangeUpperBound = 0.0
    cell(ic)%charRangeUpperBound = ' '
    select case ( node_id(root) )
    case ( n_set_one )        ! field is "/fieldname"
      son = subtree(1,root)
      call get_string ( sub_rosa(son), cell(ic)%keyword, myupcase )
      call get_string ( sub_rosa(son), cell(ic)%charValue, myupcase )
      cell(ic)%type = log_value
      cell(ic)%units = phyq_dimensionless
      cell(ic)%realValue = 1.0
    case ( n_equal, n_asg )   ! field is "name=expr"
      son = subtree(1,root)
      call get_string ( sub_rosa(son), cell(ic)%keyword, myupcase )
      do i = 2, nsons(root)
        if ( i > 2 ) then
          ic = ic + 1
          if ( ic > size(cell) ) then
            call announce_error ( root, msgnum )
            call output ( 'Number of them = ' )
            call output ( ic, advance='yes' )
            stop
          end if
          cell(ic)%keyword = ' '
          cell(ic)%charvalue = ' '
          cell(ic)%units = 0
          cell(ic)%rangeUpperBound = 0.0
          cell(ic)%charRangeUpperBound = ' '
        end if
        cell(ic)%more = nsons(root) - i
        cell(ic)%realValue = 0.0
        son = subtree(i,root)
        select case ( node_id(son) )
        case ( n_colon ) ! A range
          call expr ( son, units, value, type )
          if ( type /= str_range ) then
            cell(ic)%units = units(1)
            cell(ic)%realValue = value(1)
            cell(ic)%rangeUpperBound = value(2)
            cell(ic)%type = range
          else
            gson1 = subtree(1,son); gson2 = subtree(2,son)
            call get_string ( sub_rosa(gson1), cell(ic)%charValue, myupcase )
            call get_string ( sub_rosa(gson2), cell(ic)%charRangeUpperBound, &
                              myupcase )
            cell(ic)%type = str_range
          end if
        case ( n_dot ) ! A reference of the form X.Y
          call get_string ( sub_rosa(subtree(1,son)), cell(ic)%charValue, &
                            myupcase )
          call get_string ( sub_rosa(subtree(2,son)), &
                            cell(ic)%charRangeUpperBound, myupcase )
          cell(ic)%type = label
          cell(ic)%units = phyq_dimensionless
        case ( n_string )
          call get_string ( sub_rosa(son), cell(ic)%charValue, myupcase )
          cell(ic)%type = str_value
          cell(ic)%units = phyq_dimensionless
        case default
          call fill_cell_expr ( son, cell(ic) )
        end select
      end do
    case default              ! field is "expr"
      cell(ic)%keyword = ' '
      call fill_cell_expr ( root, cell(ic) )
    end select
    if ( toggle(gen) ) call trace_end ( 'FILL_CELL' )
  end subroutine FILL_CELL
! -----------------------------------------------  FILL_CELL_EXPR  -----
  subroutine FILL_CELL_EXPR ( ROOT, CELL )
  ! Fill a cell's RHS from expr
    integer, intent(in) :: ROOT              ! Root of subtree
    type(MlscfCell_T), intent(inout) :: CELL ! Cell to fill
    integer :: UNITS(2)                      ! Output of expr
    double precision :: VALUE(2)             ! Output of expr
    call expr ( root, units, value, cell%type )
    if ( cell%type == label ) then
      call get_string ( sub_rosa(root), cell%charValue, myupcase )
    else if ( node_id(root) == n_identifier .and. &
         ( cell%type == undeclared .or. cell%type == empty .or. &
           cell%type == enum_value) ) then
      call get_string ( sub_rosa(root), cell%charValue, myupcase )
      cell%type = enum_value
    else if ( cell%type == str_value ) then
      call get_string ( sub_rosa(root), cell%charValue, myupcase )
      cell%type = str_value
    else
      cell%units = units(1)
      cell%realValue = value(1)
    end if
  end subroutine FILL_CELL_EXPR
! -------------------------------------------------------  ONE_CF  -----
  subroutine ONE_CF ( ROOT )
  ! Analyze one configuration, with abstract syntax tree rooted at ROOT.
  ! The root should be an N_CF node (but we don't check).
    integer, intent(in) :: ROOT

    integer :: GSON1, GSON2        ! Grandsons of ROOT being analyzed
    integer :: I, N                ! Loop inductor, number of sons
    integer :: SON1, SONN          ! Sons of ROOT being analyzed
    integer :: STRING1             ! Sub_Rosa of first sons

    if ( toggle(gen) ) call trace_begin ( 'ONE_CF', root )
    isection = isection + 1
    idef = 0
    ispec = 0
    n = nsons(root)
    son1 = subtree(1,root)
    string1 = sub_rosa(son1)
    call get_string ( string1, l2cf_data%sections(isection)%MlscfSectionName, &
                      myupcase )
    do i = 2, n-1
      sonn = subtree(i,root)
      select case ( node_id(sonn) )
      case ( n_equal )           ! Spec of the form x = expr
        call fill_cell ( sonn, l2cf_data%sections(isection)%cells, idef, &
          too_many_defs )
      case ( n_named )
        ispec = ispec + 1
        if ( ispec > size(l2cf_data%sections(isection)%entries) ) then
          call announce_error ( sonn, too_many_specs )
          call output ( 'ISPEC = ' ); call output ( ispec, advance='yes' )
          stop
        end if
        gson1 = subtree(1,sonn)  ! Label
        call get_string ( sub_rosa(gson1), &
          l2cf_data%sections(isection)%entries(ispec)%MlscfLabelName, myupcase )
        gson2 = subtree(2,sonn)  ! spec_args tree
        call spec_args ( gson2 )
      case ( n_spec_args )       ! Spec of the form x,(y[=expr])+
        ispec = ispec + 1
        if ( ispec > size(l2cf_data%sections(isection)%entries) ) then
          call announce_error ( sonn, too_many_specs )
          call output ( 'ISPEC = ' ); call output ( ispec, advance='yes' )
          stop
        end if
        l2cf_data%sections(isection)%entries(ispec)%MlscfLabelName = ' '
        call spec_args ( sonn )
      case default
        call announce_error ( sonn, no_code_for )
      end select
    end do
    l2cf_data%sections(isection)%NoSectionDefs = idef
    l2cf_data%sections(isection)%NoSectionEntries = ispec
    if ( toggle(gen) ) call trace_end ( 'ONE_CF' )
  end subroutine ONE_CF
! ----------------------------------------------------  SPEC_ARGS  -----
  subroutine SPEC_ARGS ( ROOT )
  ! Analyze a specification of the form NAME (, EXPR ( = EXPR)? )+
  ! starting at ROOT
    integer, intent(in) :: ROOT

    integer :: I              ! Loop inductor
    integer :: SON            ! I'th son of "root"

    if ( toggle(gen) ) call trace_begin ( 'SPEC_ARGS', root )
    son = subtree(1,root)
    call get_string ( sub_rosa(son), &
      l2cf_data%sections(isection)%entries(ispec)%MlscfEntryName, myupcase )
    icell = 0
    do i = 2, nsons(root)
      son = subtree(i,root)
      call fill_cell ( son, &
        l2cf_data%sections(isection)%entries(ispec)%cells, icell, &
        too_many_cells )
    end do
    l2cf_data%sections(isection)%entries(ispec)%MlscfEntryNoKeys = icell
    if ( toggle(gen) ) call trace_end ( 'SPEC_ARGS' )
  end subroutine SPEC_ARGS

  end subroutine GENERATE_TABLE
end module TABLE_GENERATOR

! $Log$
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
! Revision 1.2  2000/09/01 21:42:15  vsnyder
! Add "more" field to MLSCF_CELL type
!
