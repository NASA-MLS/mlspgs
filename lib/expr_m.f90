! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module EXPR_M

! Evaluate an expression from its tree.

  implicit NONE
  private
  public :: EXPR, GetIndexFlagsFromList

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================
  ! -------------------------------------------------------  EXPR  -----
  recursive subroutine EXPR ( ROOT, UNITS, VALUE, TYPE, SCALE )
  ! Analyze an expression, return its type, units and value.

    use DECLARATION_TABLE, only: DECLARED, DECLS, EMPTY, ENUM_VALUE, &
                                 EXPRN, FUNCTION, GET_DECL, NAMED_VALUE, &
                                 NUM_VALUE, RANGE, STR_RANGE, STR_VALUE, &
                                 UNDECLARED, UNITS_NAME
    use Functions, only: F_Exp, F_Log, F_Sqrt
    use INTRINSIC, only: PHYQ_DIMENSIONLESS, PHYQ_INVALID
    use STRING_TABLE, only: FLOAT_VALUE
    use TOGGLES, only: CON, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: NODE_ID, NSONS, SUB_ROSA, SUBTREE
    use TREE_TYPES ! Everything, especially everything beginning with N_

    integer, intent(in) :: ROOT         ! Root of expression subtree
    integer, intent(out) :: UNITS(2)    ! Units of expression value -- UNITS(2)
                                        ! is PHYQ_INVALID if ROOT is not a
                                        ! range (:) operator.
    double precision, intent(out) :: VALUE(2) ! Expression value, if any
    integer, intent(out), optional :: TYPE    ! Expression type
    double precision, optional, intent(out) :: SCALE(2) ! Scale for units

    type(decls) :: DECL            ! Declaration record for "root"
    integer :: I
    integer :: ME                  ! node_id(root)
    integer :: Son1                ! First son of "root"
    integer :: STRING              ! Sub_rosa(root)
    integer :: TYPE2               ! Type of son of "root"
    integer :: UNITS2(2)           ! Units of an expression
    double precision :: VALUE2(2)  ! Value of an expression
    double precision :: SCALE2(2)  ! Units scale

    ! Error codes
    integer, parameter :: BadNode = 1 ! Unsupported tree node in expr
    integer, parameter :: NonNumeric = badNode + 1 ! Non numeric arg for func
    integer, parameter :: NotFunc = nonNumeric + 1 ! Not a function
    integer, parameter :: NotUnitless = notFunc + 1 ! Not a function
    integer, parameter :: OutOfRange = notUnitless + 1
    integer, parameter :: UnsupportedFunc = outOfRange + 1
    integer, parameter :: WrongNumArgs = unsupportedFunc + 1

    if ( toggle(con) ) call trace_begin ( 'EXPR', root )
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
    case ( n_and, n_or )
      if ( present(type) ) type = exprn
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
            call announceError ( subtree(2,root), notUnitless )
          else
            select case ( decl%units ) ! the function index in this case
            case ( f_exp )
              if ( value(1) > log(huge(value(1))) ) then
                call announceError ( subtree(2,root), outOfRange )
              else
                value(1) = exp(value(1))
              end if
            case ( f_log )
              if ( value(1) <= 0.0 ) then
                call announceError ( subtree(2,root), outOfRange )
              else
                value(1) = log(value(1))
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
        if ( type2 == num_value ) then
          value = value2(1) ** value
        else
          value = value2 ** value
        end if
      end do
    case default
      call expr ( subtree(1,root), units, value, type, scale )
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
          call expr ( subtree(2,root), units2, value2, type, scale2 )
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
            if ( me == n_plus ) then
              value = value + value2
            else !  me == n_minus
              value = value - value2
            end if
          else if ( me == n_minus ) then
            value = - value
          end if
        case ( n_mult )
          value = value * value2
          where ( units == phyq_dimensionless ) units = units2
        case ( n_div )
          value = value / value2
        case ( n_into )
          value = value2 / value
          units = units2
        case default
          call announceError ( root, badNode )
        end select
      end if
    end select
    if ( toggle(con) ) call trace_end ( 'EXPR' )
  contains
    subroutine AnnounceError ( where, what )
      use MORETREE, only: StartErrorMessage
      use OUTPUT_M, only: OUTPUT
      use String_Table, only: Display_String
      use TREE, only: DUMP_TREE_NODE
      integer, intent(in) :: Where ! Tree index
      integer, intent(in) :: What  ! Error index
      call startErrorMessage ( where )
      select case ( what )
      case ( badNode )
        call dump_tree_node ( where, 0 )
        call output ( ' is not supported.', advance='yes' )
      case ( nonNumeric )
        call output ( 'Argument of ' )
        call display_string ( string )
        call output ( ' is not numeric.', advance='yes' )
      case ( notFunc )
        call display_string ( string )
        call output ( ' is not a valid function.', advance='yes' )
      case ( notUnitless )
        call output ( 'Argument of ' )
        call display_string ( string )
        call output ( ' is not unitless.', advance='yes' )
      case ( outOfRange )
        call output ( 'Argument of ' )
        call display_string ( string )
        call output ( ' is out of range.', advance='yes' )
      case ( unsupportedFunc )
        call output ( 'Function ' )
        call display_string ( string )
        call output ( ' is not supported.', advance='yes' )
      case ( wrongNumArgs )
        call output ( 'Incorrect number of arguments for ' )
        call display_string ( string, advance='yes' )
      end select
      ! There's no way to return an error, so return something
      if ( present(type) ) type = empty
      units = PHYQ_INVALID
      value = 0.0
    end subroutine AnnounceError
  end subroutine EXPR

  ! --------------------------------------  GetIndexFlagsFromList  -----

  subroutine GetIndexFlagsFromList ( root, flags, status, lower, noError )
    ! Given the root of a numeric/numeric range array
    ! Set the flags array appropriately
    use Declaration_table, only: NUM_VALUE
    use Intrinsic, only: PHYQ_DIMENSIONLESS
    use MLSCommon, only: R8
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

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module EXPR_M

! $Log$
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
