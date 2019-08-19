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

  implicit none
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

    use Call_Stack_M, only: Stack_Depth
    use Declaration_Table, only: Allocate_Test, Deallocate_Test, Declared, &
      & Decls, Dump_Values, Empty, Enum_Value, Function, Label, Log_Value, &
      & Named_Value, Num_Value, Get_Decl, Range, Str_Range, Str_Value, &
      & Units_Name, Variable
    use Functions, only: F_Difference, F_Exp, F_Ln, F_Log, F_Log10, &
      & F_Intersection, F_Sqrt, F_Union, F_Without
    use Intrinsic, only: Data_Type_Indices, Lit_Indices, Phyq_Dimensionless, &
      & Phyq_Indices, Phyq_Invalid, L_False, L_True, T_Boolean
    use Output_M, only: Blanks, Newline, Output
    use StarterrorMessage_M, only: StarterrorMessage
    use String_Table, only: Display_String, Float_Value
    use Toggles, only: Con, Levels, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Decoration, Node_Id, Nsons, Sub_Rosa, Subtree, Tree_Node_Name
    use Tree_Types ! Everything, Especially Everything Beginning With N

    integer, intent(in) :: ROOT         ! Root of expression subtree
    integer, intent(out) :: UNITS(2)    ! Units of expression value if type ==
                                        ! num_value -- UNITS(2) is
                                        ! PHYQ_INVALID if ROOT is not a range
                                        ! (:) operator.  Type of expression
                                        ! value if type == enum_value
    double precision, intent(out) :: VALUE(2) ! Expression value, if any.  If
                                        ! TYPE == Log_Value, zero means false.
                                        ! Lit index if type is enum_value.
    integer, intent(out), optional :: TYPE    ! Expression type
    double precision, optional, intent(out) :: SCALE(2) ! Scale for units
    type(value_t), allocatable, intent(out), optional :: Values(:)

    type(decls) :: DECL            ! Declaration record for "root"
    integer :: EnumType            ! Enum data type index if MyType == enum_type
                                   ! or MyType == variable and VarType ==
                                   ! enum_type
    integer :: I, J
    integer :: ME                  ! node_id(root)
    integer :: MyType              ! Expression type as in Declaration_Table
    integer :: N                   ! Size for allocating
    logical :: OK                  ! No errors reported
    integer :: Son1                ! First son of "root"
    integer :: STRING              ! Sub_rosa(root)
    integer :: Trace = -1          ! String index for trace
    integer :: TYPE1, TYPE2        ! Type of son of "root"
    integer :: UNITS2(2)           ! Units of an expression
    double precision :: VALUE2(2)  ! Value of an expression
    integer :: VarType             ! Variable type if MyType == variable
    double precision :: SCALE2(2)  ! Units scale
    type(value_t), allocatable :: Values1(:), Values2(:), Values3(:)

    ! Error codes
    integer, parameter :: BadNode = 1 ! Unsupported tree node in expr
    integer, parameter :: DifferentShapes = badNode + 1 ! Shapes of operands differ
    integer, parameter :: DifferentTypes = differentShapes + 1 ! Types of operands differ
    integer, parameter :: NoArray = differentTypes + 1 ! No VALUES argument
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
        decl = decls(0.0d0, empty, phyq_invalid, 0, 0, null() )
      end if
      myType = decl%type
      if ( myType /= variable ) &
        & call allocate_test ( values1, 1, 'Values1', moduleName )
      select case ( myType )
      case ( enum_value )
        enumType = decoration(subtree(1,decl%tree)) ! data type index
        units = enumType
        value = decl%units                          ! lit index
                         ! What      Type     Value Units Decor
        values1 = value_t(enum_value,enumType,value,0,    0)
      case ( variable )
        call allocate_test ( values1, size(decl%values), 'Values1', moduleName )
        values1 = decl%values
        if ( decl%units == enum_value ) then
          varType = enum_value
          myType = varType
          enumType = decl%values(1)%type
          units = enumType
          value = decl%values(1)%value ! It's in values1
        else
          myType = decl%values(1)%type
          varType = decl%values(1)%units(1)
          units = varType
          value = decl%values(1)%value ! It's in values1
        end if
      case default
        units = decl%units
        value = decl%value
                         ! What  Type   value units decor
        values1 = value_t(myType,myType,value,units,0)
      end select
      if ( present(values) ) then
        ! Make sure values is allocated so we don't need to test for it
        ! after calling EXPR recursively.
        call allocate_test ( values, size(values1), 'Values', moduleName )
      else if ( decl%type == variable ) then
        if ( size(decl%values) > 1 ) call announceError ( root, noArray )
      end if
    case ( n_number ) ! ----------------------------------------------------
      string = sub_rosa(root)
      units = phyq_dimensionless
      value = float_value(string)
      myType = num_value
    case ( n_string ) ! ----------------------------------------------------
      units = phyq_dimensionless
      value = sub_rosa(root)
      myType = str_value
    case ( n_func_ref ) ! --------------------------------------------------
      son1 = subtree(1,root)
      ! Look up the function name
      string = sub_rosa(son1)
      decl = get_decl(string,function)
      if ( decl%type /= function ) call announceError ( son1, notFunc )
      if ( nsons(root) < 2 ) call announceError ( root, wrongNumArgs )
      if ( OK ) then
        call expr ( subtree(2,root), units, value, type1, values=values1 )
        myType = type1
        if ( myType == enum_value ) enumType = units(1)
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
            end if
          end if
        case ( f_difference )
          if ( nsons(root) /= 3 ) then
            call announceError ( root, wrongNumArgs )
          else
            call allocate_test ( values3, 0, 'Values3', moduleName )
            call expr ( subtree(3,root), units2, value2, type2, values=values2 )
            myType = type2
            call difference ( values1, values2 )
            call move_value ( values3 )
          end if
        case ( f_intersection )
          if ( nsons(root) /= 3 ) then
            call announceError ( root, wrongNumArgs )
          else
            call allocate_test ( values3, 0, 'Values3', moduleName )
            call expr ( subtree(3,root), units2, value2, type2, values=values2 )
            myType = type2
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
          myType = type2
          call move_value ( values3 )
        case ( f_without )
          if ( nsons(root) /= 3 ) then
            call announceError ( root, wrongNumArgs )
          else
            call allocate_test ( values3, 0, 'Values3', moduleName )
            call expr ( subtree(3,root), units2, value2, type2, values=values2 )
            myType = type2
            call without ( values1, values2 )
            call move_value ( values3 )
          end if
        case default
          call announceError ( son1, unsupportedFunc )
        end select
      end if
    case ( n_pow ) ! -------------------------------------------------------
      call expr ( subtree(nsons(root),root), units, value, type1, scale, values1 )
      myType = num_value
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
        call expr ( subtree(2,root), units, value, myType, scale, values1 )
      else
        call expr ( subtree(3,root), units, value, myType, scale, values1 )
      end if
      if ( myType == enum_value ) enumType = units(1)
    case ( n_array ) ! -------------------------------------------------
      if ( .not. present(values) ) then
        if ( nsons(root) > 1 ) then
          call announceError ( root, noArray )
        else
          call expr ( subtree(i,root), units, value, myType, scale )
          if ( myType == enum_value ) enumType = units(1)
        end if
      else
        call expr ( subtree(1,root), units, value2, type1, scale, values1 )
        myType = type1
        if ( myType == enum_value ) enumType = units(1)
        do i = 2, nsons(root)
          call expr ( subtree(i,root), units, value2, type2, scale, values2 )
          call allocate_test ( values3, size(values1) + size(values2), &
            & 'Values3', moduleName )
          values3 = [ values1, values2 ]
          call deallocate_test ( values1, 'Values1', moduleName )
          call deallocate_test ( values2, 'Values2', moduleName )
          call move_alloc ( values3, values1 )
        end do
        call set_values ( type2 ) ! Makes values1 deallocated
      end if
    case default
      call expr ( subtree(1,root), units, value, type1, scale, values1 )
      myType = type1
      if ( myType == enum_value ) enumType = units(1)
      if ( me == n_unit ) then
        decl = get_decl(sub_rosa(subtree(2,root)), units_name)
        units = decl%units
        values1%units(1) = units(1)
        values1%units(2) = units(2)
        if ( decl%value > 0.0 ) then
          value(1) = value(1) * decl%value
        else
          value(1) = value(1) - decl%value
        end if
        values1(1)%value(1) = value(1)
        if ( present(scale) ) scale(1) = decl%value
      else
        if ( nsons(root) > 1 ) &
          & call expr ( subtree(2,root), units2, value2, type2, scale2, values2 )
        select case ( me )
        case ( n_colon, n_colon_less, n_less_colon, n_less_colon_less ) ! --
          if ( size(values1) /= size(values2) ) &
            & call announceError ( root, differentShapes )
          units(2) = units2(1); value(2) = value2(1)
          values1%units(2) = units(2)
          if ( present(scale) ) scale(2) = scale2(1)
          if ( myType == num_value ) myType = range
          if ( myType == str_value ) myType = str_range
          values1%what = myType
          values1%decor = me
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
        case ( n_equal_equal, n_not_equal )
          myType = log_value
          if ( type1 /= type2 ) &
            & call announceError ( root, differentTypes )
          if ( type1 == num_value ) then
            if ( .not. check_units ( values1, values2 ) ) &
              &  call announceError ( root, wrongUnits )
          end if
          if ( OK ) then
            do j = 1, 2
              select case ( me )
              case ( n_equal_equal ) ! -------------------------------------
                values1%value(j) = merge(l_true,l_false,values1%value(j) == values2%value(j))
              case ( n_not_equal ) ! ---------------------------------------
                values1%value(j) = merge(l_true,l_false,values1%value(j) /= values2%value(j))
              end select
            end do
          end if
        case ( n_less, n_less_eq, n_greater, n_greater_eq )
          myType = log_value
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
                values1%value(j) = merge(l_true,l_false, &
                  & values1%value(j) < values2%value(j))
              case ( n_less_eq ) ! -----------------------------------------
                values1%value(j) = merge(l_true,l_false, &
                  & values1%value(j) <= values2%value(j))
              case ( n_greater ) ! -----------------------------------------
                values1%value(j) = merge(l_true,l_false, &
                  & values1%value(j) > values2%value(j))
              case ( n_greater_eq ) ! --------------------------------------
                values1%value(j) = merge(l_true,l_false, &
                  & values1%value(j) >= values2%value(j))
              end select
            end do
          end if
        case ( n_and, n_or ) ! ---------------------------------------------
          myType = log_value
          if ( .not. is_boolean(values1(1)) ) &
             & call announceError ( subtree(1,root), notLogical )
          if ( .not. is_boolean(values2(1)) ) &
             &  call announceError ( subtree(2,root), notLogical )
          if ( OK ) then
            if ( me == n_and ) then
              values1%value(1) = merge(l_true,l_false, &
                & values1%value(1)==l_true .and. values2%value(1)==l_true)
            else if ( me == n_or ) then
              values1%value(1) = merge(l_true,l_false, &
                & values1%value(1)==l_true .or. values2%value(1)==l_true)
            end if
          end if
        case ( n_not ) ! ---------------------------------------------------
          myType = log_value
          if ( .not. is_boolean(values1(1)) ) then
            call announceError ( subtree(1,root), notLogical )
          else
            values1%value(1) = merge(l_true,l_false,values1%value(1)==l_false)
          end if
        case default
          call announceError ( root, badNode )
        end select
      end if
    end select

    if ( myType == log_value ) then
      myType = enum_value
      enumType = t_boolean
      values1%what = enum_value
      values1%type = t_boolean
      value = values1%value(1)
    end if

! This doesn't work.  It demonstrates that I ought to write an expr that
! always puts its results in values, and a wrapper that takes the old
! [value,type,units] from that.
!    if ( allocated(values1) ) value = values1(1)%value

    n = 1
    if ( present(values) ) then
      if ( allocated(values) ) then
        n = size(values)
        if ( n == 0 ) then
          call deallocate_test ( values, 'Values', moduleName )
        else if ( allocated(values1) ) then
          if ( n /= size(values1) ) then
            call deallocate_test ( values, 'Values', moduleName )
          else
            values = values1
          end if
        end if
      end if
      if ( .not. allocated(values) ) then
        ! Make sure values is allocated so we don't need to test for it
        ! after calling EXPR recursively.
        if ( allocated(values1) ) then
          call allocate_test ( values, size(values1), 'Values', moduleName )
          values = values1
        else
          call allocate_test ( values, 1, 'Values', moduleName )
          select case ( myType )
          case ( variable )
                            ! What    Type    Value Units Decor
            values = value_t(variable,varType,units,0,    0)
          case ( enum_value )
                            ! What      Type     Value Units Decor
            values = value_t(enum_value,enumType,value,0,    0)
          case default
                            ! What  Type   Value Units Decor
            values = value_t(myType,myType,value,units,0)
          end select
        end if
      end if
      if ( present(type) ) then
        type = values(1)%type
        if ( type == empty ) type = myType
        if ( values(1)%what == enum_value ) then
          type = enum_value
          units = values(1)%type
        end if
      end if
      if ( allocated(values1) ) then
        values = values1
      else if ( me /= n_identifier .and. me /= n_array ) then
        select case ( values(1)%what )
        case ( enum_value, label )
          values%type = type1
        case default
          values%what = values%type
        end select
      end if
      value = values(1)%value
      units = values(1)%units
    else if ( present(type) ) then
      type = myType
    end if
    if ( myType == enum_value ) units = enumType

    if ( toggle(con) .and. levels(con) > 1 ) then
      call Blanks( Stack_Depth()+1, FillChar='_' )
      if ( me == n_identifier ) &
        & call display_string ( string, before=' Identifier ' )
      select case ( myType )
      case ( enum_value )
        call display_string ( data_type_indices(enumType), before=' enum Type = ' )
        if ( me /= n_array ) &
          & call display_string ( lit_indices(nint(value(1))), before=' Value = ' )
      case ( variable )
        if ( varType == enum_value ) then
          call display_string ( data_type_indices(enumType), &
            & before=' variable enum Type = ' )
          if ( size(decl%values) == 1 ) then
            call display_string ( lit_indices(nint(decl%values(1)%value(1))), &
              & before=' Value = ' )
          end if
        else
          call display_string ( data_type_indices(varType), &
            & before=' variable Type = ' ) 
          if ( size(decl%values) == 1 ) &
            & call display_value ( decl%values(1)%value )
        end if
      case default
        call display_string ( data_type_indices(myType), before=' Type = ' )
        call display_value ( value )
      end select
      call newLine
      if ( present(values) ) then
        if ( size(values) > 1 ) call dump_values ( values )
      end if
    end if
    if ( myType /= enum_value ) then
      call trace_end ( 'EXPR', index=n, string='Type=', &
        & stringIndex=data_type_indices(myType), cond=toggle(con) )
    else
      call trace_end ( 'EXPR', index=n, string='Type=', &
        & stringIndex=data_type_indices(enumType), cond=toggle(con) )
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
      use Tree, only: Dump_Tree_Node, Dump_Tree_Node_Name, Node_Kind, &
        & Pseudo, Sub_Rosa
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
      case ( differentTypes )
        call dump_tree_node ( where, 0 )
        call output ( ' has operands with different types.', advance='yes' )
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
        call dump_tree_node_name ( where, before='Argument of ' )
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

    subroutine Display_Units ( Index )
      ! Display the units and scale of a numeric value
      integer, intent(in) :: Index
      if ( units(index) /= 0 ) then
        call display_string ( phyq_indices(units(index)), before=' ' )
      else
        call output ( units(index), before=' unit ' )
      end if
      if ( present(scale) ) call output ( scale(index), before=' Scale = ' )
    end subroutine Display_Units

    subroutine Display_Value ( Value )
      ! Display the value as a number, a number range, a string, or a
      ! string range
      double precision, intent(in) :: Value(2)
      call output ( ' Value = ' )
      select case ( myType )
      case ( num_value )
        call output ( value(1) )
        call display_units ( 1 )
      case ( range )
        call output ( value(1) )
        call display_units ( 1 )
        call output ( value(2), before=' : ' )
        call display_units ( 2 )
      case ( str_range )
        call display_string ( nint(value(1)), strip=.false. )
      case ( str_value )
        call display_string ( nint(value(1)), strip=.false. )
        call display_string ( nint(value(2)), strip=.false., before=' : ' )
      end select
    end subroutine Display_Value

    subroutine Intersection ( V1, V2 )
      ! Add elements that are both V1 and V2 to Values3
      ! Add elements of V2 that are not in V1 to Values3
      type(value_t), intent(in) :: V1(:), V2(:)
      integer :: I
      do i = 1, size(v1)
        if ( any(v1(i) == v2) ) call add_to_set ( v1(i) )
      end do
    end subroutine Intersection

    logical function Is_Boolean ( Value )
      type(value_t), intent(in) :: Value
      is_boolean = value%what == enum_value .and. value%type == t_boolean
    end function Is_Boolean

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
      if ( allocated(values1) ) then
        if ( present(values) ) then
          call move_alloc ( values1, values )
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
    use Intrinsic, only: Phyq_Invalid
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
    use Intrinsic, only: Phyq_Dimensionless
    use MLSKinds, only: R8
    use Tree, only: Node_ID, Subtree, Nsons
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
! Revision 2.36  2019/08/19 22:02:55  pwagner
! Light housekeeping
!
! Revision 2.35  2014/03/20 01:40:16  vsnyder
! Unified types in Intrinsic, repair some tracing
!
! Revision 2.34  2014/03/05 01:06:50  vsnyder
! Repair blunder that caused a^b^c...x^y^z to be evaluated as a^b^c...x^y^a.
!
! Revision 2.33  2014/02/27 02:27:40  vsnyder
! Corrections to units and ranges
!
! Revision 2.32  2014/02/21 19:26:03  vsnyder
! Extensive work for variables and enumeration-type results
!
! Revision 2.31  2014/01/11 01:41:02  vsnyder
! Decruftification
!
! Revision 2.30  2014/01/08 21:11:22  vsnyder
! More type checking.  Better handling of arrays.  Allow == and /= for
! types other than numeric.
!
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
