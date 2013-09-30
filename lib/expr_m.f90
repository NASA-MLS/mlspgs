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

  implicit NONE
  private
  public :: EXPR, EXPR_CHECK, GetIndexFlagsFromList

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================
  ! -------------------------------------------------------  EXPR  -----
  recursive subroutine EXPR ( ROOT, UNITS, VALUE, TYPE, SCALE )
  ! Analyze an expression, return its type, units and value.

    use DECLARATION_TABLE, only: DECLARED, DECLS, EMPTY, ENUM_VALUE, &
                                 FUNCTION, GET_DECL, LOG_VALUE, &
                                 NAMED_VALUE, NUM_VALUE, RANGE, STR_RANGE, &
                                 STR_VALUE, UNDECLARED, UNITS_NAME
    use Functions, only: F_Exp, F_Ln, F_Log, F_Log10, F_Sqrt
    use INTRINSIC, only: PHYQ_DIMENSIONLESS, PHYQ_INVALID
    use StartErrorMessage_m, only: StartErrorMessage
    use STRING_TABLE, only: FLOAT_VALUE
    use TOGGLES, only: CON, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: NODE_ID, NSONS, SUB_ROSA, SUBTREE
    use TREE_TYPES ! Everything, especially everything beginning with N_

    integer, intent(in) :: ROOT         ! Root of expression subtree
    integer, intent(out) :: UNITS(2)    ! Units of expression value -- UNITS(2)
                                        ! is PHYQ_INVALID if ROOT is not a
                                        ! range (:) operator.
    double precision, intent(out) :: VALUE(2) ! Expression value, if any.  If
                                        ! TYPE == Log_Value, zero means false.
    integer, intent(out), optional :: TYPE    ! Expression type
    double precision, optional, intent(out) :: SCALE(2) ! Scale for units

    type(decls) :: DECL            ! Declaration record for "root"
    integer :: I
    integer :: ME                  ! node_id(root)
    integer :: Son1                ! First son of "root"
    integer :: STRING              ! Sub_rosa(root)
    integer :: Trace = -1          ! String index for trace
    integer :: TYPE1, TYPE2        ! Type of son of "root"
    integer :: UNITS2(2)           ! Units of an expression
    double precision :: VALUE2(2)  ! Value of an expression
    double precision :: SCALE2(2)  ! Units scale

    ! Error codes
    integer, parameter :: BadNode = 1 ! Unsupported tree node in expr
    integer, parameter :: NonNumeric = badNode + 1 ! Non numeric arg for func
    integer, parameter :: NotFunc = nonNumeric + 1 ! Not a function
    integer, parameter :: NotLogical = notFunc + 1 ! Not logical
    integer, parameter :: NotUnitless = notLogical + 1 ! Not unitless
    integer, parameter :: NotUnitlessArg = notUnitless + 1 ! Not unitless
    integer, parameter :: OutOfRange = notUnitlessArg + 1
    integer, parameter :: UnsupportedFunc = outOfRange + 1
    integer, parameter :: WrongNumArgs = unsupportedFunc + 1
    integer, parameter :: wrongUnits = wrongNumArgs + 1

    call trace_begin ( trace, 'EXPR', root, cond=toggle(con) )
    units = (/ phyq_dimensionless, phyq_invalid /)     ! default
    value = 0.0d0                                      ! default
    if ( present(scale) ) scale = 1.0d0                ! default
    me = node_id(root)
    select case ( me )
    case ( n_identifier )
      string = sub_rosa(root)
      if ( declared(string) ) then
        decl = get_decl(string, named_value)
        if ( decl%type == empty ) decl = get_decl(string, enum_value)
      else
        decl = decls(0.0d0, undeclared, phyq_invalid, 0, 0 )
      end if
      if ( present(type) ) type = decl%type
      units = decl%units
      value = decl%value
    case ( n_number )
      string = sub_rosa(root)
      units = phyq_dimensionless
      value = float_value(string)
      if ( present(type) ) type = num_value
    case ( n_string )
      if ( present(type) ) type = str_value
      units = phyq_dimensionless
      value = sub_rosa(root)
    case ( n_func_ref )
      son1 = subtree(1,root)
      ! Look up the function name
      string = sub_rosa(son1)
      decl = get_decl(string,function)
      if ( decl%type /= function ) then
        call announceError ( son1, notFunc )
      else
        if ( nsons(root) /= 2 ) then
          call announceError ( root, wrongNumArgs )
        else
          call expr ( subtree(2,root), units, value, type2 )
          if ( type2 /= num_value ) then
            call announceError ( subtree(2,root), nonNumeric )
          else if ( units(1) /= phyq_dimensionless ) then
            !??? Does the tree checker already check this?
            call announceError ( subtree(2,root), notUnitlessArg )
          else
            select case ( decl%units ) ! the function index in this case
            case ( f_exp )
              if ( value(1) > log(huge(value(1))) ) then
                call announceError ( subtree(2,root), outOfRange )
              else
                value(1) = exp(value(1))
              end if
            case ( f_ln, f_log )
              if ( value(1) <= 0.0 ) then
                call announceError ( subtree(2,root), outOfRange )
              else
                value(1) = log(value(1))
              end if
            case ( f_log10 )
              if ( value(1) <= 0.0 ) then
                call announceError ( subtree(2,root), outOfRange )
              else
                value(1) = log10(value(1))
              end if
            case ( f_sqrt )
              if ( value(1) < 0.0 ) then
                call announceError ( subtree(2,root), outOfRange )
              else
                value(1) = sqrt(value(1))
              end if
            case default
              call announceError ( son1, unsupportedFunc )
            end select
            if ( present(type) ) type = num_value
          end if
        end if
      end if
    case ( n_pow )
      value = 1.0
      do i = nsons(root), 1, -1 ! Power operator is right associative
        call expr ( subtree(i,root), units, value2, type2 )
        if ( any(units /= phyq_dimensionless ) ) then
            call announceError ( subtree(2,root), notUnitless )
        else if ( type2 == num_value ) then
          value = value2(1) ** value
        else
          value = value2 ** value
        end if
      end do
    case default
      call expr ( subtree(1,root), units, value, type1, scale )
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
          & call expr ( subtree(2,root), units2, value2, type2, scale2 )
        select case ( me )
        case ( n_colon, n_colon_less, n_less_colon, n_less_colon_less )
          units(2) = units2(1); value(2) = value2(1)
          if ( present(scale) ) scale(2) = scale2(1)
          if ( present(type) ) then
            if ( type == num_value ) type = range
            if ( type == str_value ) type = str_range
          end if
        case ( n_plus, n_minus )
          if ( nsons(root) > 1 ) then
            if ( any(units /= units2) ) call announceError ( root, wrongUnits )
            if ( me == n_plus ) then
              value = value + value2
            else !  me == n_minus
              value = value - value2
            end if
          else if ( me == n_minus ) then
            value = - value
          end if
        case ( n_mult )
          if ( .not. any(units == phyq_dimensionless) .and. &
             & .not. any(units2 == phyq_dimensionless) ) &
             & call announceError ( root, wrongUnits )
          value = value * value2
          where ( units == phyq_dimensionless ) units = units2
        case ( n_div )
          if ( .not. any(units == phyq_dimensionless) .and. &
             & .not. any(units2 == phyq_dimensionless) ) &
             & call announceError ( root, wrongUnits )
          value = value / value2
        case ( n_into )
          if ( .not. any(units == phyq_dimensionless) .and. &
             & .not. any(units2 == phyq_dimensionless) ) &
             & call announceError ( root, wrongUnits )
          value = value2 / value
          units = units2
        case ( n_less, n_less_eq, n_greater, n_greater_eq, n_equal_equal, n_not_equal )
          if ( present(type) ) type = log_value
          if ( type1 /= num_value ) then
            call announceError ( subtree(1,root), nonNumeric )
            if ( type2 /= num_value ) &
              call announceError ( subtree(2,root), nonNumeric )
          else if ( type2 /= num_value ) then
            call announceError ( subtree(2,root), nonNumeric )
          else if ( any(units /= units2) ) then
            call announceError ( root, wrongUnits )
          else
            select case ( me )
            case ( n_less )
              value(1) = merge(1.0,0.0,value(1) < value2(1) )
            case ( n_less_eq )
              value(1) = merge(1.0,0.0,value(1) <= value2(1) )
            case ( n_greater )
              value(1) = merge(1.0,0.0,value(1) > value2(1) )
            case ( n_greater_eq )
              value(1) = merge(1.0,0.0,value(1) >= value2(1) )
            case ( n_equal_equal )
              value(1) = merge(1.0,0.0,value(1) == value2(1) )
            case ( n_not_equal )
              value(1) = merge(1.0,0.0,value(1) /= value2(1) )
            end select
          end if
        case ( n_and, n_or )
          if ( present(type) ) type = log_value
          if ( type1 /= log_value ) then
            call announceError ( subtree(1,root), notLogical )
            if ( type2 /= log_value ) &
              call announceError ( subtree(2,root), notLogical )
          else if ( type2 /= log_value ) then
            call announceError ( subtree(2,root), notLogical )
          else if ( me == n_and ) then
            value(1) = value(1) * value2(1)
          else if ( me == n_or ) then
            value(1) = max(value(1), value2(1))
          end if
        case ( n_not )
          if ( present(type) ) type = log_value
          if ( type1 /= log_value ) then
            call announceError ( subtree(1,root), notLogical )
          else
            value(1) = 1 - nint(value(1))
          end if
        case default
          call announceError ( root, badNode )
        end select
      end if
    end select
    call trace_end ( 'EXPR', index=type, cond=toggle(con) )
  contains
    subroutine AnnounceError ( where, what )
      use Output_m, only: Output
      use String_Table, only: Display_String
      use Tree, only: Dump_Tree_Node
      integer, intent(in) :: Where ! Tree index
      integer, intent(in) :: What  ! Error index
      call startErrorMessage ( where )
      select case ( what )
      case ( badNode )
        call dump_tree_node ( where, 0 )
        call output ( ' is not supported.', advance='yes' )
      case ( nonNumeric )
        call display_string ( string, before='Argument of ' )
        call output ( ' is not numeric.', advance='yes' )
      case ( notFunc )
        call display_string ( string )
        call output ( ' is not a valid function.', advance='yes' )
      case ( notLogical )
        call display_string ( string, before='Argument of ' )
        call output ( ' is not logical.', advance='yes' )
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
