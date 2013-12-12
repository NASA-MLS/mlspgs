! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module EXPR_M

! Evaluate an expression from its tree.

  use Declaration_table, only: Value_T ! Needed for Value_ procedures
                                       ! which cannot be internal until
                                       ! more compilers allow generic interfaces
                                       ! for internal procedures.

  implicit NONE
  private
  public :: EXPR, EXPR_CHECK, GetIndexFlagsFromList

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  include "Value_T_Interfaces.f9h"

contains ! ====     Public Procedures     ==============================
  ! -------------------------------------------------------  EXPR  -----
  recursive subroutine EXPR ( ROOT, UNITS, VALUE, TYPE, SCALE, VALUES  )
  ! Analyze an expression, return its type, units and value.
  ! If its result is an array, return the values in VALUES.

    use DECLARATION_TABLE, only: Allocate_Test, Deallocate_Test, DECLARED, &
      DECLS, Dump_Values, EMPTY, ENUM_VALUE, FUNCTION, LABEL, LOG_VALUE, &
      NAMED_VALUE, NUM_VALUE, GET_DECL, RANGE, STR_RANGE, STR_VALUE, &
      Type_Names, Type_Name_Indices, UNDECLARED, UNITS_NAME, Variable
    use Functions, only: F_Difference, F_Exp, F_Ln, F_Log, F_Log10, &
      F_Intersection, F_Sqrt, F_Union, F_Without
    use INTRINSIC, only: PHYQ_DIMENSIONLESS, PHYQ_Indices, PHYQ_INVALID, T_Unknown
    use Output_m, only: NewLine, Output
    use StartErrorMessage_m, only: StartErrorMessage
    use STRING_TABLE, only: Display_String, FLOAT_VALUE
    use TOGGLES, only: CON, LEVELS, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: Decoration, NODE_ID, NSONS, SUB_ROSA, SUBTREE, Tree_Node_Name
    use TREE_TYPES ! Everything, especially everything beginning with N_

    integer, intent(in) :: ROOT         ! Root of expression subtree
    integer, intent(out) :: UNITS(2)    ! Units of expression value -- UNITS(2)
                                        ! is PHYQ_INVALID if ROOT is not a
                                        ! range (:) operator.
    double precision, intent(out) :: VALUE(2) ! Expression value, if any.  If
                                        ! TYPE == Log_Value, zero means false.
    integer, intent(out), optional :: TYPE    ! Expression type
    double precision, optional, intent(out) :: SCALE(2) ! Scale for units
    type(value_t), allocatable, intent(out), optional :: Values(:)

    type(decls) :: DECL            ! Declaration record for "root"
    integer :: I, J
    integer :: ME                  ! node_id(root)
    integer :: N                   ! Size for allocating
    logical :: OK                  ! No errors reported
    integer :: Son1                ! First son of "root"
    integer :: STRING              ! Sub_rosa(root)
    integer :: Trace = -1          ! String index for trace
    integer :: TYPE1, TYPE2        ! Type of son of "root"
    integer :: UNITS2(2)           ! Units of an expression
    double precision :: VALUE2(2)  ! Value of an expression
    double precision :: SCALE2(2)  ! Units scale
    type(value_t), allocatable :: Values1(:), Values2(:), Values3(:)

    ! Error codes
    integer, parameter :: BadNode = 1 ! Unsupported tree node in expr
    integer, parameter :: DifferentShapes = badNode + 1 ! Shapes of operands differ
    integer, parameter :: NoArray = differentShapes + 1 ! No VALUES argument
    integer, parameter :: NonNumeric = noArray + 1 ! Non numeric arg for func
    integer, parameter :: NotFunc = nonNumeric + 1 ! Not a function
    integer, parameter :: NotLogical = notFunc + 1 ! Not logical
    integer, parameter :: NotScalar = notLogical + 1  ! Not scalar
    integer, parameter :: NotUnitless = notScalar + 1 ! Not unitless
    integer, parameter :: NotUnitlessArg = notUnitless + 1 ! Not unitless
    integer, parameter :: OutOfRange = notUnitlessArg + 1
    integer, parameter :: UnsupportedFunc = outOfRange + 1
    integer, parameter :: WrongNumArgs = unsupportedFunc + 1
    integer, parameter :: wrongUnits = wrongNumArgs + 1

    call trace_begin ( trace, 'EXPR', root, string=tree_node_name(root), &
      & cond=toggle(con) )
    units = (/ phyq_dimensionless, phyq_invalid /)     ! default
    value = 0.0d0                                      ! default
    if ( present(scale) ) scale = 1.0d0                ! default
    if ( present(values) ) call deallocate_test ( values, 'Values', moduleName )
    OK = .true.
    me = node_id(root)
    select case ( me )
    case ( n_identifier ) ! ------------------------------------------------
      string = sub_rosa(root)
      if ( declared(string) ) then
        decl = get_decl ( string, [ named_value, variable, enum_value, label ] )
      else
        decl = decls(0.0d0, undeclared, phyq_invalid, 0, 0, null() )
      end if
      if ( present(type) ) type = decl%type
      units = decl%units
      value = decl%value
      if ( present(values) ) then
        if ( decl%type == variable .and. allocated(decl%values) ) then
          n = size(decl%values)
          call allocate_test ( values, size(decl%values), 'Values', moduleName )
          values = decl%values
        else
          call allocate_test ( values, 1, 'Values', moduleName )
          values = value_t(decl%type,decl%type,value,units)
          select case ( decl%type )
          case ( enum_value, label )
            ! Make sure values is allocated so we don't need to test for it
            ! after calling EXPR recursively.
            type1 = decl%tree ! dt_def node if enum_value, first son is type name
                              ! spec_args node if label, first son is spec name,
                              ! decoration of first son is spec_def node
            if ( decl%type == label ) type1 = decoration(subtree(1,type1))
          end select
        end if
      end if
    case ( n_number ) ! ----------------------------------------------------
      string = sub_rosa(root)
      units = phyq_dimensionless
      value = float_value(string)
      if ( present(type) ) type = num_value
    case ( n_string ) ! ----------------------------------------------------
      if ( present(type) ) type = str_value
      units = phyq_dimensionless
      value = sub_rosa(root)
    case ( n_func_ref ) ! --------------------------------------------------
      son1 = subtree(1,root)
      ! Look up the function name
      string = sub_rosa(son1)
      decl = get_decl(string,function)
      if ( decl%type /= function ) call announceError ( son1, notFunc )
      if ( nsons(root) < 2 ) call announceError ( root, wrongNumArgs )
      if ( OK ) then
        call expr ( subtree(2,root), units, value, type1, values=values1 )
        if ( present(type) ) type = type1
        select case ( decl%units ) ! the function index in this case
        case ( f_exp, f_ln, f_log, f_log10, f_sqrt ) ! One-argument numeric
          if ( nsons(root) /= 2 ) call announceError ( root, wrongNumArgs )
          if ( OK ) then
            n = 0
            if ( type1 == num_value ) n = 1
            if ( type1 == range ) n = 2
            if ( n == 0 ) then
              call announceError ( subtree(2,root), nonNumeric )
            else if ( .not. check_units(values1, unit=phyq_dimensionless) ) then
              !??? Does the tree checker already check this?
              call announceError ( subtree(2,root), notUnitlessArg )
            else
              select case ( decl%units ) ! the function index in this case
              case ( f_exp )
                if ( any(values1%value(1) > log(huge(values1%value(1)))) ) then
                  call announceError ( subtree(2,root), outOfRange )
                else
                  values1%value(1) = exp(values1%value(1))
                end if
              case ( f_ln, f_log )
                if ( any(values1%value(1) <= 0.0) ) then
                  call announceError ( subtree(2,root), outOfRange )
                else
                  values1%value(1) = log(values1%value(1))
                end if
              case ( f_log10 )
                if ( any(values1%value(1) <= 0.0) ) then
                  call announceError ( subtree(2,root), outOfRange )
                else
                  values1%value(1) = log10(values1%value(1))
                end if
              case ( f_sqrt )
                if ( any(values1%value(1) < 0.0) ) then
                  call announceError ( subtree(2,root), outOfRange )
                else
                  values1%value(1) = sqrt(values1%value(1))
                end if
              end select
              call set_values ( num_value )
              if ( present(type) ) type = type1
            end if
          end if
        case ( f_difference )
          if ( nsons(root) /= 3 ) then
            call announceError ( root, wrongNumArgs )
          else
            call allocate_test ( values3, 0, 'Values3', moduleName )
            call expr ( subtree(3,root), units2, value2, type2, values=values2 )
            call difference ( values1, values2 )
            call move_value ( values3 )
          end if
        case ( f_intersection )
          if ( nsons(root) /= 3 ) then
            call announceError ( root, wrongNumArgs )
          else
            call allocate_test ( values3, 0, 'Values3', moduleName )
            call expr ( subtree(3,root), units2, value2, type2, values=values2 )
            call intersection ( values1, values2 )
            call move_value ( values3 )
          end if
        case ( f_union )
          call allocate_test ( values3, 0, 'Values3', moduleName )
          do j = 1, size(values1)
            call add_to_set ( values1(j) )
          end do
          do i = 3, nsons(root)
            call expr ( subtree(i,root), units2, value2, type2, values=values2 )
            do j = 1, size(values2)
              call add_to_set ( values2(j) )
            end do
          end do
          call move_value ( values3 )
        case ( f_without )
          if ( nsons(root) /= 3 ) then
            call announceError ( root, wrongNumArgs )
          else
            call allocate_test ( values3, 0, 'Values3', moduleName )
            call expr ( subtree(3,root), units2, value2, type2, values=values2 )
            call without ( values1, values2 )
            call move_value ( values3 )
          end if
        case default
          call announceError ( son1, unsupportedFunc )
        end select
      end if
    case ( n_pow ) ! -------------------------------------------------------
      call expr ( subtree(1,root), units, value, type1, scale, values1 )
      do i = nsons(root)-1, 1, -1 ! Power operator is right associative
        call expr ( subtree(i,root), units, value2, type2, scale, values2 )
        if ( size(values1) /= size(values2) ) then
          call announceError ( root, differentShapes )
        else
          if ( .not. check_units(values2, unit=phyq_dimensionless) ) then
            call announceError ( subtree(i,root), notUnitless )
          else
            do j = 1, size(values2)
              if ( values2(j)%type == num_value ) then
                values1(j)%value = values2(j)%value(1) ** values1(j)%value
              else
                values1(j)%value = values2(j)%value ** values1(j)%value
              end if
            end do
          end if
          call set_values ( type2 )
        end if
      end do
    case ( n_cond ) ! ------------------------------------------------------
      ! expr ? expr ! expr
      ! First subtree is known to be logical
      call expr ( subtree(1,root), units, value, type1, scale, values1 )
      if ( size(values1) /= 1 ) &
        & call announceError ( subtree(1,root), notScalar )
      if ( value(1) /= 0 ) then ! subtree 1 is true
        call expr ( subtree(2,root), units, value, type, scale, values1 )
      else
        call expr ( subtree(3,root), units, value, type, scale, values1 )
      end if
    case default
      call expr ( subtree(1,root), units, value, type1, scale, values1 )
      if ( present(type) ) type = type1
      if ( me == n_unit ) then
        decl = get_decl(sub_rosa(subtree(2,root)), units_name)
        units = decl%units
        if ( decl%value > 0.0 ) then
          value(1) = value(1) * decl%value
        else
          value(1) = value(1) - decl%value
        end if
        if ( present(scale) ) scale(1) = decl%value
      else
        if ( nsons(root) > 1 ) &
          & call expr ( subtree(2,root), units2, value2, type2, scale2, values2 )
        select case ( me )
        case ( n_colon, n_colon_less, n_less_colon, n_less_colon_less ) ! --
          if ( size(values1) /= size(values2) ) &
            & call announceError ( root, differentShapes )
          units(2) = units2(1); value(2) = value2(1)
          if ( present(scale) ) scale(2) = scale2(1)
          if ( present(type) ) then
            if ( type == num_value ) type = range
            if ( type == str_value ) type = str_range
          end if
        case ( n_plus, n_minus ) ! -----------------------------------------
          if ( nsons(root) > 1 ) then
            if ( .not. check_units ( values1, values2 ) ) &
               & call announceError ( root, wrongUnits )
            if ( me == n_plus ) then
              values1 = values1 + values2
            else !  me == n_minus
              values1 = values1 - values2
            end if
            call set_values ( type1, type2 )
          else if ( me == n_minus ) then
            values1 = -values1
            call set_values ( type1 )
          end if
        case ( n_mult ) ! --------------------------------------------------
          ! At least one operand must be entirely dimensionless.
          ! We should do this on an element-by-element bases.
          if ( .not. check_units(values1, unit=phyq_dimensionless) .and. &
             & .not. check_units(values2, unit=phyq_dimensionless) ) &
             & call announceError ( root, wrongUnits )
          values1 = values1 * values2
          do j = 1, 2
            where ( values1%units(j) == phyq_dimensionless ) &
              values1%units(j) = values2%units(j)
          end do
        case ( n_div ) ! ---------------------------------------------------
          if ( .not. check_units(values2, unit=phyq_dimensionless) ) &
             & call announceError ( root, wrongUnits )
          values1 = values1 / values2
        case ( n_into ) ! --------------------------------------------------
          if ( .not. check_units(values1, unit=phyq_dimensionless) ) &
             & call announceError ( root, wrongUnits )
          values1 = values2 / values1
          do j = 1, 2
            values1%units(j) = values2%units(j)
          end do
        case ( n_less, n_less_eq, n_greater, n_greater_eq, n_equal_equal, n_not_equal )
          if ( present(type) ) type = log_value
          if ( type1 /= num_value ) &
             & call announceError ( subtree(1,root), nonNumeric )
          if ( type2 /= num_value ) &
             & call announceError ( subtree(2,root), nonNumeric )
          if ( .not. check_units ( values1, values2 ) ) &
             &  call announceError ( root, wrongUnits )
          if ( OK ) then
            do j = 1, 2
              select case ( me )
              case ( n_less ) ! --------------------------------------------
                values1%value(j) = merge(1.0,0.0,values1%value(j) < values2%value(j))
              case ( n_less_eq ) ! -----------------------------------------
                values1%value(j) = merge(1.0,0.0,values1%value(j) <= values2%value(j))
              case ( n_greater ) ! -----------------------------------------
                values1%value(j) = merge(1.0,0.0,values1%value(j) > values2%value(j))
              case ( n_greater_eq ) ! --------------------------------------
                values1%value(j) = merge(1.0,0.0,values1%value(j) >= values2%value(j))
              case ( n_equal_equal ) ! -------------------------------------
                values1%value(j) = merge(1.0,0.0,values1%value(j) == values2%value(j))
              case ( n_not_equal ) ! ---------------------------------------
                values1%value(j) = merge(1.0,0.0,values1%value(j) /= values2%value(j))
              end select
            end do
            call set_values ( log_value )
          end if
        case ( n_and, n_or ) ! ---------------------------------------------
          if ( present(type) ) type = log_value
          if ( type1 /= log_value ) &
             & call announceError ( subtree(1,root), notLogical )
          if ( type2 /= log_value ) &
             &  call announceError ( subtree(2,root), notLogical )
          if ( OK ) then
            if ( me == n_and ) then
              values1%value(1) = values1%value(1) * values2%value(1)
            else if ( me == n_or ) then
              values1%value(1) = max(values1%value(1), values2%value(1))
            end if
            call set_values ( log_value )
          end if
        case ( n_not ) ! ---------------------------------------------------
          if ( present(type) ) type = log_value
          if ( type1 /= log_value ) then
            call announceError ( subtree(1,root), notLogical )
          else
            values1%value(1) = 1 - nint(values1%value(1))
            call set_values ( log_value )
          end if
        case ( n_array ) ! -------------------------------------------------
          if ( .not. present(values) ) then
            if ( nsons(root) > 1 ) then
              call announceError ( root, noArray )
            else
              call expr ( subtree(i,root), units, value, type, scale )
            end if
          else
            call deallocate_test ( values3, 'Values3', moduleName )
            call allocate_test ( values1, 1, 'Values1', moduleName )
            call expr ( subtree(1,root), units, value2, type2, scale, values2 )
            values1(1) = value_t(empty,type,value,units)
            call set_values ( type2 ) ! Makes values1 deallocated
            do i = 2, nsons(root)
              call expr ( subtree(i,root), units, value2, type2, scale, values2 )
              n = size(values) + size(values2)
              call allocate_test ( values1, n, 'Values1', moduleName )
              values1(1:size(values)) = values
              values1(size(values)+1:n) = values2
              call set_values ( type2 ) ! Makes values1 deallocated
            end do
          end if
        case default
          call announceError ( root, badNode )
        end select
      end if
    end select

    n = 1
    if ( present(values) ) then
      if ( allocated(values) ) then
        n = size(values)
        if ( n == 0 ) &
          call deallocate_test ( values, 'Values', moduleName )
      end if
      if ( .not. allocated(values) ) then
        ! Make sure values is allocated so we don't need to test for it
        ! after calling EXPR recursively.
        call allocate_test ( values, 1, 'Values', moduleName )
        values = value_t(empty,type,value,units)
      end if
      if ( present(type) ) type = values(1)%type
      select case ( values(1)%type )
      case ( enum_value, label )
        values%type = type1
      case default
        values%what = values%type
      end select
      value = values(1)%value
      units = values(1)%units
    end if

    if ( toggle(con) .and. levels(con) > 1 ) then
      if ( node_id(root) == n_identifier ) then
        call display_string ( string, before='Value = ' )
      else
        call output ( value(1), before='Value = ' )
        if ( units(1) /= 0 ) then
          call display_string ( phyq_indices(units(1)), before=' ' )
        else
          call output ( units(1), before=' unit ' )
        end if
        call output ( value(2), before=' ' )
        if ( units(2) /= 0 ) then
          call display_string ( phyq_indices(units(2)), before=' ' )
        else
          call output ( units(2), before=' unit ' )
        end if
      end if
      if ( present(type) ) &
        & call display_string ( type_name_indices(type), before=' Type = ' )
      if ( present(scale) ) then
        call output ( scale(1), before=' Scale = ' )
        call output ( scale(2), before=' ' )
      end if
      call newLine
      if ( present(values) ) then
        if ( size(values) > 1 ) call dump_values ( values )
      end if
    end if
    if ( present(type) ) then
      call trace_end ( 'EXPR', index=n, string='Type=' //trim(type_names(type)), &
        & cond=toggle(con) )
    else
      call trace_end ( 'EXPR', index=n, cond=toggle(con) )
    end if

  contains

    subroutine Add_To_Set ( Value )
      ! Add Value to Values3, which is assumed to be allocated
      type(value_t), intent(in) :: Value
      type(value_t), allocatable :: Temp(:)

      if ( any(value == values3) ) return ! No duplicates
      call allocate_test ( temp, size(values3)+1, 'Temp', moduleName )
      temp(:size(values3)) = values3
      call deallocate_test ( values3, 'Values3', moduleName )
      call move_alloc ( temp, values3 )
      values3(size(values3)) = value
    end subroutine Add_To_Set

    subroutine AnnounceError ( where, what )
      use Output_m, only: Output
      use String_Table, only: Display_String
      use Tree, only: Dump_Tree_Node, Dump_Tree_Node_Name, Internal, &
        & Node_Kind, Pseudo, Sub_Rosa
      integer, intent(in) :: Where ! Tree index
      integer, intent(in) :: What  ! Error index
      integer :: NodeIs
      OK = .false.
      call startErrorMessage ( where )
      select case ( what )
      case ( badNode )
        call dump_tree_node ( where, 0 )
        call output ( ' is not supported.', advance='yes' )
      case ( differentShapes )
        call dump_tree_node ( where, 0 )
        call output ( ' has operands with different shapes.', advance='yes' )
      case ( noArray )
        call dump_tree_node ( where, 0 )
        call output ( ' EXPR invoked for array but VALUES argument not present.', &
          advance='yes' )
      case ( nonNumeric )
        nodeIs = node_kind ( where )
        if ( nodeIs == pseudo ) then
          call display_string ( sub_rosa(where) )
        else
          call dump_tree_node_name ( where, before='Argument of ' )
        end if
        call output ( ' is not numeric.', advance='yes' )
      case ( notFunc )
        call display_string ( string )
        call output ( ' is not a valid function.', advance='yes' )
      case ( notLogical )
        call display_string ( string, before='Argument of ' )
        call output ( ' is not logical.', advance='yes' )
      case ( notScalar )
        call dump_tree_node ( where, 0 )
        call output ( ' does not have a scalar value.', advance='yes' )
      case ( notUnitless )
        call output ( 'Operands of ' )
        call dump_tree_node ( where, 0 )
        call output ( ' are not unitless.', advance='yes' )
      case ( notUnitlessArg )
        call display_string ( string, before='Argument of ' )
        call output ( ' is not unitless.', advance='yes' )
      case ( outOfRange )
        call display_string ( string, before='Argument of ' )
        call output ( ' is out of range.', advance='yes' )
      case ( unsupportedFunc )
        call output ( 'Function ' )
        call display_string ( string )
        call output ( ' is not supported.', advance='yes' )
      case ( wrongNumArgs )
        call output ( 'Incorrect number of arguments for ' )
        call display_string ( string, advance='yes' )
      case ( wrongUnits )
        call output ( 'Units of operands of ' )
        call dump_tree_node ( where, 0 )
        call output ( ' are not compatible.', advance='yes' )
      end select
      ! There's no way to return an error, so return something
      if ( present(type) ) type = empty
      units = PHYQ_INVALID
      value = 0.0
    end subroutine AnnounceError

    logical function Check_Units ( V1, V2, Unit ) result ( R )
      ! If Unit is present return true if all units in V1 are equal
      ! to Unit, and if V2 is also present its units are also equal
      ! to Unit.
      ! If Unit is not present, assume both V1 and V2 are present, and
      ! return true if all their units are equal.
      type(value_t), intent(in) :: V1(:)
      type(value_t), intent(in), optional :: V2(:)
      integer, intent(in), optional :: Unit
      integer :: J
      r = .true.
      if ( present(unit) ) then
        if ( present(v2) ) then
          do j = 1, 2
            r = r .and. all(v1%units(j) == unit) .and. all(v2%units(j) == unit)
          end do
        else
          do j = 1, 2
            r = r .and. all(v1%units(j) == unit)
          end do
        end if
      else
        do j = 1, 2
          r = r .and. all(v1%units(j) == v2%units(j))
        end do
      end if
    end function Check_Units

    subroutine Difference ( V1, V2 )
      ! Add elements of V1 that are not in V2 to Values3
      ! Add elements of V2 that are not in V1 to Values3
      type(value_t), intent(in) :: V1(:), V2(:)
      integer :: I

      do i = 1, size(v1)
        if ( all(v1(i) /= v2) ) call add_to_set ( v1(i) )
      end do
      do i = 1, size(v2)
        if ( all(v1 /= v2(i)) ) call add_to_set ( v2(i) )
      end do
    end subroutine Difference

    subroutine Intersection ( V1, V2 )
      ! Add elements that are both V1 and V2 to Values3
      ! Add elements of V2 that are not in V1 to Values3
      type(value_t), intent(in) :: V1(:), V2(:)
      integer :: I
      do i = 1, size(v1)
        if ( any(v1(i) == v2) ) call add_to_set ( v1(i) )
      end do
    end subroutine Intersection

    subroutine Move_Value ( Value )
      ! Move Value to Values
      type(value_t), intent(inout), allocatable :: Value(:)
      call move_alloc ( value, values )
    end subroutine Move_Value

    subroutine Set_Values ( Type1, Type2 )
      ! Set the Values argument if it is present and Values1 is allocated
      ! with size other than 1.  If Values1 is allocated with size 1 put
      ! Values1%type in type and Values1%value in value.  Deallocate
      ! Values1 and Values2.
      integer, intent(in) :: Type1
      integer, intent(in), optional :: Type2
      integer :: Type
      type = type1
      if ( present(type2) ) type = the_type ( type1, type2 )
      if ( allocated(values1) ) then
        if ( present(values) ) then
          call move_alloc ( values1, values )
          values%type = type
          values%what = decl%type
          select case ( decl%type )
          case ( enum_value )
            values%type = decl%units
          case ( variable )
            values%type = decl%units
            if ( allocated(decl%values) ) units = decl%values(1)%units(1)
          end select
        else
          value = values1(1)%value
          call deallocate_test ( values1, 'Values1', moduleName )
        end if
      end if
      call deallocate_test ( values2, 'Values2', moduleName )
    end subroutine Set_Values

    integer function The_Type ( Type1, Type2 )
      integer, intent(in) :: Type1, Type2
      if ( type1 == type2 ) then
        the_type = type1
      else if ( type1 == num_value ) then
        the_type = range
      else
        the_type = str_range
      end if
    end function The_Type

    subroutine Without ( V1, V2 )
      ! Add elements of V1 that are not in V2 to Values3
      type(value_t), intent(in) :: V1(:), V2(:)
      integer :: I

      do i = 1, size(v1)
        if ( all(v1(i) /= v2) ) call add_to_set ( v1(i) )
      end do
    end subroutine Without

  end subroutine EXPR

  ! -------------------------------------------------  EXPR_CHECK  -----
  subroutine EXPR_CHECK ( ROOT, UNITS, VALUE, NEED, ERROR, TYPE, SCALE )
  ! Analyze an expression, return its type, units and value.  Check that
  ! its units are one of the units in NEED.  ERROR = true if not.
    use INTRINSIC, only: PHYQ_INVALID
    integer, intent(in) :: ROOT         ! Root of expression subtree
    integer, intent(out) :: UNITS(2)    ! Units of expression value -- UNITS(2)
                                        ! is PHYQ_INVALID if ROOT is not a
                                        ! range (:) operator.
    double precision, intent(out) :: VALUE(2) ! Expression value, if any
    integer, intent(in) :: NEED(:)      ! Needed units
    logical, intent(out) :: ERROR       ! "Wrong units"
    integer, intent(out), optional :: TYPE    ! Expression type
    double precision, optional, intent(out) :: SCALE(2) ! Scale for units
    integer :: I

    call expr ( root, units, value, type, scale )
    error = .true.
    do i = 1, size(need)
      if ( units(1) == need(i) ) then
        error = .false.
        exit
      end if
    end do
    if ( units(2) /= phyq_invalid .and. .not. error ) then
      error = .true.
      do i = 1, size(need)
        if ( units(2) == need(i) ) then
          error = .false.
          exit
        end if
      end do
    end if

  end subroutine EXPR_CHECK

  ! --------------------------------------  GetIndexFlagsFromList  -----

  subroutine GetIndexFlagsFromList ( root, flags, status, lower, noError )
    ! Given the root of a numeric/numeric range array
    ! Set the flags array appropriately
    use Declaration_table, only: NUM_VALUE
    use Intrinsic, only: PHYQ_DIMENSIONLESS
    use MLSKinds, only: R8
    use Tree, only: Node_ID, Subtree, nsons
    use Tree_Types, only: N_colon_less, N_less_colon, N_less_colon_less
    integer, intent(in) :: ROOT         ! Tree node
    logical, dimension(:), intent(inout) :: FLAGS ! Result
    integer, intent(out) :: STATUS      ! Error flag, 0=success
    integer, intent(in), optional :: LOWER ! Lower bound for result
    logical, intent(in), optional :: NOERROR ! If set don't give bounds errors

    ! Local variables
    integer :: I,J                      ! Loop counters
    real(r8), dimension(2) :: VALUE     ! From expr
    integer, dimension(2) :: UNITS      ! From expr
    integer :: TYPE                     ! From expr
    integer :: RANGE_LOW, RANGE_HI      ! Range for flags
    logical :: MYNOERROR                ! Copy of noError
    integer :: MYLOWER                  ! Copy of lower
    integer :: SON                      ! Tree node

    ! Executable code
    flags = .false.
    status = 0
    myNoError = .false.
    if ( present ( noError ) ) myNoError = noError
    mylower = 1
    if ( present ( lower ) ) myLower = lower

    do j = 2, nsons(root)
      son = subtree ( j, root )
      call expr ( son, units, value, type )
      do i = 1, merge(1,2,type==num_value)
        if ( units(i) /= phyq_dimensionless ) then
          status = 1
          return
        end if
      end do
      range_low = nint(value(1))
      range_hi = nint(value(merge(1,2,type==num_value)))
      select case ( node_id(son) )
      case ( n_colon_less )
        range_hi = range_hi - 1
      case ( n_less_colon )
        range_low = range_low + 1
      case ( n_less_colon_less )
        range_low = range_low + 1
        range_hi = range_hi - 1
      end select
      if ( .not. myNoError .and. &
        & ( range_low < myLower .or. range_hi > myLower+size(flags)-1 ) ) then
        status = 1
        return
      end if
      range_low = max ( range_low, myLower )
      range_hi = min ( range_hi, myLower+size(flags)-1 )
      flags ( range_low-myLower+1 : range_hi-myLower+1 ) = .true.
    end do
  end subroutine GetIndexFlagsFromList

  include "Value_T_Implementations.f9h"

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module EXPR_M

! $Log$
! Revision 2.29  2013/12/12 02:03:41  vsnyder
! Add ability to substitute variable values
!
! Revision 2.28  2013/10/16 01:12:02  vsnyder
! Include Null() in for allocatable component in Decls() constructor.
! Don't try to dump Values if it's not present.
!
! Revision 2.27  2013/10/15 23:55:20  pwagner
! NAG demanded this change to intent in Move_Value
!
! Revision 2.26  2013/10/11 01:47:21  vsnyder
! Revision for array results seems to work
!
! Revision 2.25  2013/09/30 23:59:24  vsnyder
! Move StartErrorMessage from include to module
!
! Revision 2.24  2013/09/21 00:35:50  vsnyder
! More trace output
!
! Revision 2.23  2013/09/19 23:27:38  vsnyder
! More careful units checking
!
! Revision 2.22  2013/08/30 16:44:09  pwagner
! Trying to fix bug in call trace usage
!
! Revision 2.21  2013/08/30 03:56:02  vsnyder
! Revise use of trace_begin and trace_end
!
! Revision 2.20  2012/05/07 23:00:57  vsnyder
! StartErrorMessage moved to include to avoid a circular dependence
! between expr_m and MoreTree
!
! Revision 2.19  2012/05/05 00:11:51  vsnyder
! Add support for 'not' operator
!
! Revision 2.18  2012/04/24 20:38:47  vsnyder
! Get kinds from MLSKinds instead of MLSCommon
!
! Revision 2.17  2012/03/12 18:36:11  vsnyder
! Add ln as synonym for log, add log10
!
! Revision 2.16  2011/04/19 01:59:43  vsnyder
! Support == and /= relational operators too
!
! Revision 2.15  2011/04/18 19:33:26  vsnyder
! Add support for relational operators and boolean-valued expressions
!
! Revision 2.14  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.13  2008/09/04 20:03:09  vsnyder
! Add PRINT statement in not_used_here
!
! Revision 2.12  2005/08/04 02:55:02  vsnyder
! Add Expr_Check
!
! Revision 2.11  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.10  2004/05/28 23:45:09  vsnyder
! Get units from either operand of *, second operand of \\
!
! Revision 2.9  2004/05/28 23:13:46  vsnyder
! Add power (^) operator, log, exp and sqrt functions
!
! Revision 2.8  2004/05/28 00:57:25  vsnyder
! Move GetIndexFlagsFromList from MoreTree to Expr_m
!
! Revision 2.7  2004/01/17 03:04:48  vsnyder
! Provide for functions in expressions
!
! Revision 2.6  2004/01/14 18:32:58  vsnyder
! Stuff for Algebra module
!
! Revision 2.5  2002/10/08 00:09:09  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2002/10/02 00:44:08  vsnyder
! Add optional SCALE argument
!
! Revision 2.3  2001/11/27 00:54:37  vsnyder
! Implement (partially) open ranges
!
! Revision 2.2  2001/04/09 20:59:57  vsnyder
! Subtract negative scale factors instead of multiplying
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
