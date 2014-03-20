! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Vector_Qty_Expr_m

  ! Evaluate an expression for which the value is a vector quantity
  ! or a number.

  ! Descriptions of expression results.  If Dot the result is a vector
  ! quantity.  If T_Boolean or T_Numeric the result is a number, being
  ! either L_True or L_False in the T_Boolean case.  If negative, an error
  ! occurred.

  use Intrinsic, only: L_False, L_True, Dot => T_A_dot_B, T_Boolean, T_Numeric

  implicit NONE
  private

  public :: Dot, Vector_Qty_Expr

  ! Error codes
  integer, private, parameter :: Bad_Operator = 1
  integer, private, parameter :: Bad_Terminal = bad_operator + 1
  integer, private, parameter :: Incompatible = bad_terminal + 1
  integer, private, parameter :: NotFunc = incompatible + 1
  integer, private, parameter :: OutOfRange = notFunc + 1
  integer, private, parameter :: UnsupportedFunc = outOfRange + 1
  integer, private, parameter :: WrongNumArgs = unsupportedFunc + 1
  integer, private, parameter :: WrongUnits = wrongNumArgs + 1

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================

  ! --------------------------------------------  Vector_Qty_Expr  -----
  recursive integer function Vector_Qty_Expr ( Root, Vectors, Qty, Number, &
    & Units, ForWhom ) result ( Stat )

    ! Evaluate an expression.  The result is
    !   Dot => result is a vector quantity, and Qty is that quantity;
    !          this is always a new quantity, and needs to be destroyed.
    !          Number is undefined.
    !   T_Numeric => result is a number, and Number is its value.  Qty is
    !          undefined.
    !   T_Boolean => result is boolean and Number is either L_True or L_False
    !  <0 =>   an error occurred, and both Qty and Number are undefined; the
    !          value is the negative of the error code.
    ! Range (colon, etc.) and power (^) are not supported.  The only
    ! functions supported are exp, ln, log, log10, and sqrt.
    ! Boolean-valued operators set values to l_true for true, l_false for false

    use Declaration_Table, only: Decls, Function, Get_Decl, Units_Name
    use Functions, only: F_Exp, F_Ln, F_Log, F_Log10, F_Sqrt
    use Intrinsic, only: PHYQ_Dimensionless, PHYQ_Invalid
    use QuantityTemplates, only: NullifyQuantityTemplate, QuantitiesAreCompatible
    use String_Table, only: Float_Value
    use Toggles, only: Gen, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    use Tree, only: Decoration, Node_ID, NSons, Sub_Rosa, Subtree
    use Tree_Types, only: Last_Tree_Node, N_And, N_Div, N_Dot, N_Equal_Equal, &
      N_Func_Ref, N_Greater, N_Greater_Eq, N_Identifier, N_Into, N_Less,      &
      N_Less_eq, N_Minus, N_Mult, N_Not, N_Not_Equal, N_Null, N_Number, N_Or, &
      N_Plus, N_String, N_Unit
    use VectorsModule, only: CloneVectorQuantity, &
      & DestroyVectorQuantityMask, &
      & DestroyVectorQuantityValue, GetVectorQtyByTemplateIndex, &
      & Vector_t, VectorValue_t
    integer, intent(in) :: Root               ! Tree node index
    type(vector_t), intent(in) :: Vectors(:)  ! Vectors database
    type(vectorValue_t), intent(out) :: Qty   ! Result quantity
    double precision, intent(out) :: Number   ! Result number
    integer, intent(out) :: Units             ! PHYQ_... if numeric result,
                                              ! else PHYQ_Invalid
    character(*), intent(in), optional :: ForWhom

    integer, parameter :: Neither = 0, Both = 1, Left = 2, Right = 3 ! Values for Case
    integer :: Case                      ! Which operands are vector quantities
    type(decls) :: Decl                  ! Of a function reference
    integer :: Me = -1                   ! String index for trace cacheing
    integer :: Qty_Index                 ! in the vector
    type(vectorValue_t) :: Qty_1, Qty_2  ! Vector quantity operands
    integer :: Son1, Son2                ! Function name, arg tree nodes
    integer :: Stat_1, Stat_2            ! Function result statuses
    integer :: String                    ! Function name
    integer :: Types(n_null+1:last_tree_node)
    integer :: Units_1, Units_2
    double precision :: Value_1, Value_2 ! Numeric operands
    integer :: Vector_Index              ! in the Vectors database
    integer :: What                      ! Tree node ID

    data types(n_and)         / t_boolean /
    data types(n_dot)         / dot /
    data types(n_div)         / t_numeric /
    data types(n_equal_equal) / t_boolean /
    data types(n_greater_eq)  / t_boolean /
    data types(n_greater)     / t_boolean /
    data types(n_into)        / t_numeric /
    data types(n_less_eq)     / t_boolean /
    data types(n_less)        / t_boolean /
    data types(n_minus)       / t_numeric /
    data types(n_mult)        / t_numeric /
    data types(n_not_equal)   / t_boolean /
    data types(n_not)         / t_boolean /
    data types(n_number)      / t_numeric /
    data types(n_or)          / t_boolean /
    data types(n_plus)        / t_numeric /

    call trace_begin ( me, 'Vector_Qty_Expr', root, cond=toggle(gen) )
    what = node_id(root)
    if ( any(what == (/ n_and,        n_div,    n_dot,     n_equal_equal, &
                      & n_greater,    n_into,   n_less_eq, n_less,        &
                      & n_greater_eq, n_minus,  n_mult,    n_not_equal,   &
                      & n_not,        n_number, n_or,      n_plus /) )  ) &
      & stat=types(what)
    call nullifyQuantityTemplate ( qty%template )
    units = phyq_invalid

    if ( nsons(root) <= 1 ) then
      select case ( what )
      case ( n_number ) ! ----------------------------------------------
        number = float_value(sub_rosa(root))
        units = phyq_dimensionless
      case ( n_not ) ! -------------------------------------------------
        stat = vector_qty_expr ( subtree(1,root), vectors, qty, value_1, units_1 )
        select case ( stat )
        case ( dot )
          qty%values = 1 - qty%values
        case ( t_boolean )
          number = merge(l_false,l_true,value_1==l_true)
        end select
      case ( n_minus ) ! -----------------------------------------------
        stat = vector_qty_expr ( subtree(1,root), vectors, qty, value_1, units )
        select case ( stat )
        case ( dot )
          qty%values = - qty%values
        case ( t_numeric )
          number = - value_1
        end select
      case ( n_plus ) ! ------------------------------------------------
        stat = vector_qty_expr ( subtree(1,root), vectors, qty, value_1, units )
        select case ( stat )
        case ( dot )
          qty%values = qty%values
        case ( t_numeric )
          number = value_1
        end select
      case default
        stat = announce_error ( root, bad_operator )
      end select
    else
      select case ( what )
      case ( n_func_ref ) ! --------------------------------------------
        son1 = subtree(1,root)
        ! Look up the function name
        string = sub_rosa(son1)
        decl = get_decl(string,function)
        if ( decl%type /= function ) then
          stat = announce_error ( son1, notFunc )
        else
          if ( nsons(root) /= 2 ) then
            stat = announce_error ( root, wrongNumArgs )
          else
            son2 = subtree(2,root)
            stat = vector_qty_expr ( son2, vectors, qty, value_1, units )
            if ( stat /= dot .and. units /= phyq_dimensionless ) then
              stat = announce_error ( son1, wrongUnits )
    go to 9
            end if
            select case ( decl%units )
            case ( f_exp ) ! ...........................................
              select case ( stat )
              case ( dot )
                if ( any(qty%values > log(huge(qty%values(1,1)))) ) then
                  stat = announce_error ( son2, outOfRange )
                else
                  qty%values = exp(qty%values)
                end if
              case ( t_numeric )
                if ( number > log(huge(number)) ) then
                  stat = announce_error ( son2, outOfRange, number=number )
                else
                  number = exp(value_1)
                end if
              end select
            case ( f_ln, f_log ) ! .....................................
              select case ( stat )
              case ( dot )
                if ( any(qty%values <= 0.0) ) then
                  stat = announce_error ( son2, outOfRange )
                else
                  qty%values = log(qty%values)
                end if
              case ( t_numeric )
                if ( number < 0.0 ) then
                  stat = announce_error ( son2, outOfRange, number=number )
                else
                  number = log(value_1)
                end if
              end select
            case ( f_log10 ) ! .........................................
              select case ( stat )
              case ( dot )
                if ( any(qty%values <= 0.0) ) then
                  stat = announce_error ( son2, outOfRange )
                else
                  qty%values = log10(qty%values)
                end if
              case ( t_numeric )
                if ( number < 0.0 ) then
                  stat = announce_error ( son2, outOfRange, number=number )
                else
                  number = log10(value_1)
                end if
              end select
            case ( f_sqrt ) ! ..........................................
              select case ( stat )
              case ( dot )
                if ( any(qty%values < 0.0) ) then
                  stat = announce_error ( son2, outOfRange )
                else
                  qty%values = sqrt(qty%values)
                end if
              case ( t_numeric )
                if ( number < 0.0 ) then
                  stat = announce_error ( son2, outOfRange, number=number )
                else
                  number = sqrt(value_1)
                end if
              end select
            case default
              stat = announce_error ( son1, unsupportedFunc )
            end select
          end if
        end if
      case ( n_less, n_less_eq, n_greater, n_greater_eq, n_equal_equal, &
           & n_not_equal, n_and, n_or, n_plus, n_minus, n_mult, n_div,  &
           & n_into ) ! ------------------------------------------------
        stat_1 = vector_qty_expr ( subtree(1,root), vectors, qty_1, value_1, units_1 )
        units = units_1
        stat = stat_1
        stat_2 = vector_qty_expr ( subtree(2,root), vectors, qty_2, value_2, units_2 )
        if ( stat_1 == dot .and. stat_2 == dot ) then
          if ( .not. QuantitiesAreCompatible(qty_1%template, qty_2%template) ) then
            stat = announce_error ( root, incompatible )
    go to 9
          end if
  !??? Do we want to check masks here?
          call cloneVectorQuantity ( qty, qty_1, options='d' )
          qty%label = 0
          if ( qty_1%template%name /= qty_2%template%name ) qty%template%name = 0
          case = both
        else if ( stat_1 == dot ) then
          call cloneVectorQuantity ( qty, qty_1, options='d' )
          qty%label = 0
          case = left
        else if ( stat_2 == dot ) then
          call cloneVectorQuantity ( qty, qty_2, options='d' )
          qty%label = 0
          stat = dot
          case = right
        else
          case = neither
          select case ( what )
          case ( n_less, n_less_eq, n_greater, n_greater_eq, n_equal_equal, &
            & n_not_equal, n_plus, n_minus )
            if ( units_1 /= units_2 ) stat = announce_error ( root, wrongUnits )
          end select
        end if
        if ( case /= both ) qty%template%name = 0
        select case ( what )
        case ( n_less ) ! ----------------------------------------------
          select case ( case )
          case ( neither )
            number = merge(1, 0, value_1 < value_2)
          case ( both )
            where ( qty_1%values < qty_2%values )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          case ( left )
            where ( qty_1%values < value_2 )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          case ( right )
            where ( value_1 < qty_2%values )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          end select
        case ( n_less_eq ) ! -------------------------------------------
          select case ( case )
          case ( neither )
            number = merge(1, 0, value_1 <= value_2)
          case ( both )
            where ( qty_1%values <= qty_2%values )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          case ( left )
            where ( qty_1%values <= value_2 )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          case ( right )
            where ( value_1 <= qty_2%values )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          end select
        case ( n_greater ) ! -------------------------------------------
          select case ( case )
          case ( neither )
            number = merge(1, 0, value_1 > value_2)
          case ( both )
            where ( qty_1%values > qty_2%values )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          case ( left )
            where ( qty_1%values > value_2 )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          case ( right )
            where ( value_1 > qty_2%values )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          end select
        case ( n_greater_eq ) ! ----------------------------------------
          select case ( case )
          case ( neither )
            number = merge(1, 0, value_1 >= value_2)
          case ( both )
            where ( qty_1%values > qty_2%values )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          case ( left )
            where ( qty_1%values > value_2 )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          case ( right )
            where ( value_1 > qty_2%values )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          end select
        case ( n_equal_equal ) ! ---------------------------------------
          select case ( case )
          case ( neither )
            number = merge(1, 0, value_1 == value_2)
          case ( both )
            where ( qty_1%values == qty_2%values )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          case ( left )
            where ( qty_1%values == value_2 )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          case ( right )
            where ( value_1 == qty_2%values )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          end select
        case ( n_not_equal ) ! -----------------------------------------
          select case ( case )
          case ( neither )
            number = merge(1, 0, value_1 /= value_2)
          case ( both )
            where ( qty_1%values /= qty_2%values )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          case ( left )
            where ( qty_1%values /= value_2 )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          case ( right )
            where ( value_1 /= qty_2%values )
              qty%values = 1
            elsewhere
              qty%values = 0
            end where
          end select
        case ( n_and ) ! -----------------------------------------------
          select case ( case )
          case ( neither )
            number = merge(l_true,l_false,value_1==l_true .and. value_2==l_true)
          case ( both )
            qty%values = merge(l_true,l_false,qty_1%values==l_true .and. &
                                              qty_2%values==l_true)
          case ( left )
            qty%values = merge(l_true,l_false,qty_1%values==l_true .and. &
                                              value_2==l_true)
          case ( right )
            qty%values = merge(l_true,l_false,value_1==l_true .and. &
                                              qty_2%values==l_true)
          end select
        case ( n_or ) ! ------------------------------------------------
          select case ( case )
          case ( neither )
            number = merge(l_true,l_false,value_1==l_true .or. value_2==l_true)
          case ( both )
            qty%values = merge(l_true,l_false,qty_1%values==l_true .or. &
                                              qty_2%values==l_true)
          case ( left )
            qty%values = merge(l_true,l_false,qty_1%values==l_true .or. &
                                              value_2==l_true)
          case ( right )
            qty%values = merge(l_true,l_false,value_1==l_true .or. &
                                              qty_2%values==l_true)
          end select
        case ( n_plus ) ! ----------------------------------------------
          select case ( case )
          case ( neither )
            number = value_1 + value_2
          case ( both )
            qty%values = qty_1%values + qty_2%values
          case ( left )
            qty%values = qty_1%values + value_2
          case ( right )
            qty%values = value_1 + qty_2%values
          end select
        case ( n_minus ) ! ---------------------------------------------
          select case ( case )
          case ( neither )
            number = value_1 - value_2
          case ( both )
            qty%values = qty_1%values - qty_2%values
          case ( left )
            qty%values = qty_1%values - value_2
          case ( right )
            qty%values = value_1 - qty_2%values
          end select
        case ( n_mult ) ! ----------------------------------------------
          select case ( case )
          case ( neither )
            if ( units_1 /= phyq_dimensionless .and. units_2 /= phyq_dimensionless ) &
              stat = announce_error ( root, wrongUnits )
            units = merge(units_1, units_2, units_1 /= phyq_dimensionless)
            number = value_1 * value_2
          case ( both )
            qty%values = qty_1%values * qty_2%values
          case ( left )
            qty%values = qty_1%values * value_2
          case ( right )
            qty%values = value_1 * qty_2%values
          end select
        case ( n_div ) ! -----------------------------------------------
          select case ( case )
          case ( neither )
            if ( units_2 == phyq_dimensionless ) then
              units = units_1
            else if ( units_1 == units_2 ) then
              units = phyq_dimensionless
            else
              stat = announce_error ( root, wrongUnits )
            end if
            number = value_1 / value_2
          case ( both )
            qty%values = qty_1%values / qty_2%values
          case ( left )
            qty%values = qty_1%values / value_2
          case ( right )
            qty%values = value_1 / qty_2%values
          end select
        case ( n_into ) ! ----------------------------------------------
          select case ( case )
          case ( neither )
            if ( units_1 == phyq_dimensionless ) then
              units = units_2
            else if ( units_1 == units_2 ) then
              units = phyq_dimensionless
            else
              stat = announce_error ( root, wrongUnits )
            end if
            number = value_2 / value_1
          case ( both )
            qty%values = qty_2%values / qty_1%values
          case ( left )
            qty%values = value_2 / qty_1%values
          case ( right )
            qty%values = qty_2%values / value_1
          end select
        end select
        if ( stat_1 == dot ) then
          call destroyVectorQuantityMask ( qty_1 )
          call destroyVectorQuantityValue ( qty_1 )
        end if
        if ( stat_2 == dot ) then
          call destroyVectorQuantityMask ( qty_2 )
          call destroyVectorQuantityValue ( qty_2 )
        end if
      case ( n_dot ) ! -------------------------------------------------
        vector_index = decoration(decoration(subtree(1,root)))
        qty_index = decoration(decoration(decoration(subtree(2,root))))
        call nullifyQuantityTemplate ( qty%template )
        call cloneVectorQuantity ( qty, &
          & GetVectorQtyByTemplateIndex( &
          & vectors(vector_index), qty_index ), options='d' )
      case ( n_identifier, n_string ) ! --------------------------------
        stat = announce_error ( root, bad_terminal )
      case ( n_unit ) ! ------------------------------------------------
        stat = vector_qty_expr ( subtree(1,root), vectors, qty, value_1, units )
        decl = get_decl(sub_rosa(subtree(2,root)), units_name)
        units = decl%units
        if ( decl%value > 0.0 ) then
          number = value_1 * decl%value
        else
          number = value_1 - decl%value
        end if
      case default ! ---------------------------------------------------
        stat = announce_error ( root, bad_operator )
      end select

    end if

  9 call trace_end ( 'Vector_Qty_Expr', cond=toggle(gen) )

  contains

    integer function Announce_Error ( Where, What, Number )
      use Lexer_Core, only: Where_t
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error
      use String_Table, only: Get_String, String_Length
      use Tree, only: Node_id, Subtree, Sub_Rosa, Tree_Text, Where_At=>Where
      integer, intent(in) :: Where ! Tree node index
      integer, intent(in) :: What  ! Error code
      double precision, intent(in), optional :: Number

      type(where_t) :: Where_is
      integer :: L, N
      character(len=127) :: Name, QtyNames(2,2), String

      announce_error = -what
      where_is = where_at(where)
      write ( string, 1 ) where_is%source/256, mod(where_is%source,256)
      1 format ( 'At line ', i0, ', column ', i0 )
      l = len_trim(string) + 1
      if ( where_is%file /= 0 ) then
        string(l:l+3) = ' in '
        call get_string ( where_is%file, string(l+5:), strip=.true. )
        l = len_trim(string) + 1
      end if
      select case ( what )
      case ( bad_operator )
        call get_string ( tree_text(node_id(where)), name )
        n = string_length ( tree_text(node_id(where)) )
        write ( string(l:), 3 ) name(:n)
      3 format ( ', cannot handle ', a, ' in this context.' )
      case ( bad_terminal )
        call get_string ( sub_rosa(where), name )
        n = string_length ( tree_text(node_id(where)) )
        write ( string, 3 ) name(:n)
      case ( incompatible )
        call get_string ( subtree(1,subtree(1,where)), qtyNames(1,1) )
        call get_string ( subtree(2,subtree(1,where)), qtyNames(2,1) )
        call get_string ( subtree(1,subtree(2,where)), qtyNames(1,2) )
        call get_string ( subtree(2,subtree(2,where)), qtyNames(2,2) )
        write ( string(l:), 4 ) trim(qtyNames(1,1)), trim(qtyNames(2,1)), &
                              & trim(qtyNames(1,2)), trim(qtyNames(2,2))
      4 format ( ', the quantities "', a, '.', a, '" and "', a, '.', a, &
          & '" are incompatible.' )
      case ( notFunc )
        call get_string ( sub_rosa(where), name )
        write ( string(l:), 5 ) trim(name)
      5 format ( ', "', a, '" is not a function name.' )
      case ( outOfRange )
        if ( present(number) ) then
          write ( string(l:), 6 ) number
      6   format ( ', the function argument ', g15.7, ' is out of range.' )
        else
          string(l:) = ', the function argument is out of range.'
        end if
      case ( unsupportedFunc )
        call get_string ( sub_rosa(where), name )
        write ( string(l:), 8 ) trim(name)
      8 format ( ', "', a, '" is not a supported function.' )
      case ( wrongNumArgs )
        string(l:) = ', the function has the wrong number of arguments.'
      case ( wrongUnits )
        string(l:) = ', the units are incorrect.'
      end select

      if ( present(forWhom) ) then
        call MLSMessage ( MLSMSG_Error, forWhom, trim(string) )
      else
        call MLSMessage ( MLSMSG_Error, moduleName, trim(string) )
      end if

    end function Announce_Error

  end function Vector_Qty_Expr

! ====     Private Procedures     ======================================

!--------------------------- end bloc --------------------------------------
 logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Vector_Qty_Expr_m

! $Log$
! Revision 2.5  2014/03/20 01:39:47  vsnyder
! Unified types in Intrinsic
!
! Revision 2.4  2013/09/24 23:27:14  vsnyder
! Use Get_Where or Print_Source to start error messages
!
! Revision 2.3  2013/09/21 00:36:23  vsnyder
! Cannonball polishing
!
! Revision 2.2  2013/09/19 23:24:11  vsnyder
! Remove debugging print
!
! Revision 2.1  2013/09/19 23:23:10  vsnyder
! Initial commit
!
