! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module EXPR_M

! Evaluate an expression from its tree.

  use DECLARATION_TABLE, only: DECLARED, DECLS, EMPTY, ENUM_VALUE, &
                               EXPRN, GET_DECL, NAMED_VALUE, &
                               NUM_VALUE, RANGE, STR_RANGE, STR_VALUE, &
                               UNDECLARED, UNITS_NAME
  use INTRINSIC, only: PHYQ_DIMENSIONLESS, PHYQ_INVALID
  use STRING_TABLE, only: FLOAT_VALUE
  use TOGGLES, only: CON, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: NODE_ID, NSONS, SUB_ROSA, SUBTREE
  use TREE_TYPES ! Everything, especially everything beginning with N_
  implicit NONE
  private
  public :: EXPR

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================
  ! -------------------------------------------------------  EXPR  -----
  recursive subroutine EXPR ( ROOT, UNITS, VALUE, TYPE )
  ! Analyze an expression, return its type, units and value.
    integer, intent(in) :: ROOT         ! Root of expression subtree
    integer, intent(out) :: UNITS(2)    ! Units of expression value -- UNITS(2)
                                        ! is PHYQ_INVALID if ROOT is not a
                                        ! range (:) operator.
    double precision, intent(out) :: VALUE(2)! Expression value, if any
    integer, intent(out), optional :: TYPE        ! Expression type

    type(decls) :: DECL            ! Declaration record for "root"
    integer :: ME                  ! node_id(root)
    integer :: STRING              ! Sub_rosa(root)
    integer :: UNITS2(2)           ! Units of an expression
    double precision :: VALUE2(2)  ! Value of an expression

    if ( toggle(con) ) call trace_begin ( 'EXPR', root )
    units = (/ phyq_dimensionless, phyq_invalid /)     ! default
    value = 0.0d0                  ! default
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
    case default
      call expr ( subtree(1,root), units, value, type )
      if ( me == n_unit ) then
        decl = get_decl(sub_rosa(subtree(2,root)), units_name)
        units = decl%units
        if ( decl%value > 0.0 ) then
          value(1) = value(1) * decl%value
        else
          value(1) = value(1) - decl%value
        end if
      else
        if ( nsons(root) > 1 ) &
          call expr ( subtree(2,root), units2, value2, type )
        select case ( me )
        case ( n_colon )
          units(2) = units2(1); value(2) = value2(1)
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
        case ( n_div )
          value = value / value2
        case default
          ! Shouldn't get here -- presumably checked already
        end select
      end if
    end select
    if ( toggle(con) ) call trace_end ( 'EXPR' )
  end subroutine EXPR
end module EXPR_M

! $Log$
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
