! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Evaluate_Variable_m

  implicit NONE
  private
  public :: Define_Variable_As_String, Evaluate_Variable

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================

  subroutine Define_Variable_As_String ( Variable_Name, String )
    use Declaration_Table, only: Allocate_Test, Redeclare, Str_Value, Value_t, &
      & Variable
    use Symbol_Table, only: Enter_Terminal
    use Symbol_Types, only: T_Identifier, T_String
    character(*), intent(in) :: Variable_Name, String
    integer :: Name, Value
    type(value_t), allocatable :: The_Value(:)

    name = enter_terminal ( trim(adjustl(variable_name)), t_identifier )
    value = enter_terminal ( trim(adjustl(string)), t_string )
    call allocate_test ( the_value, 1, 'The_Value', moduleName )
    ! Don't use T_String for the type here; it's a token index here.
    the_value(:) = value_t(str_value,str_value,dble(value))
    !              string  value       type      units     tree
    call redeclare ( name, name+0.0d0, variable, str_value, 0, &
      & values=the_value ) ! Declares it if it's not declared
  end subroutine Define_Variable_As_String

  subroutine Evaluate_Variable ( Root, Son )

! Evaluate a variable defined by <name> := <expr>.
! Put its value in the Values component of its declaration.

    use Declaration_Table, only: Allocate_Test, Deallocate_Test, Dump_Values, &
      & Redeclare, Value_t, Variable
    use Expr_m, only: Expr
    use Output_m, only: Output
    use String_Table, only: Display_String
    use Toggles, only: Gen, Levels, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    use Tree, only: Nsons, Subtree, Sub_Rosa

    integer, intent(in) :: Root  ! The := tree node of the definition
    integer, intent(in), optional :: Son  ! Do only this son

    integer :: Details           ! Dump expression values
    integer :: First             ! First son to do, or only if present(son)
    integer :: I
    integer :: Me = -1           ! String index for trace cacheing
    integer :: Name              ! String index  of name of variable
    integer :: Type
    integer :: Units(2), Units2(2) ! From expr
    double precision :: Value(2) ! From expr
    type(value_t), allocatable :: Values1(:), Values2(:), Values3(:)

    call trace_begin ( me, 'Evaluate_Variable', root, cond=toggle(gen) )
    details = merge(levels(gen),0,toggle(gen))
    name = sub_rosa(subtree(1,root))
    first = 2
    if ( present(son) ) first = son
    call expr ( subtree(first,root), units, value, type, values=values1 )
    if ( details > 1 ) then
      call output ( first-1, before='Element ' )
      call display_string ( name, before=' of value of ', advance='yes' )
      call dump_values ( values1, details=9 )
    end if
    if ( .not. present(son) ) then
      do i = 3, nsons(root)
        call expr ( subtree(i,root), units2, value,  type, values=values2 )
        if ( details > 1 ) then
          call output ( i-1, before='Element ' )
          call display_string ( name, before=' of value of ', advance='yes' )
          call dump_values ( values2, details=9 )
        end if
        call allocate_test ( values3, size(values1)+size(values2), &
          & 'Values3', moduleName )
        values3 = [ values1, values2 ]
        call deallocate_test ( values1, 'Values1', moduleName )
        call deallocate_test ( values2, 'Values2', moduleName )
        call move_alloc ( values3, values1 )
      end do
    end if
    if ( details > 0 ) then
      call display_string ( name, before='Final value of ', advance='yes' )
      call dump_values ( values1, details=details )
    end if
    !              string  value       type     units  tree
    call redeclare ( name, name+0.0d0, variable, type, root, &
      & values=values1 ) ! Declares it if it wasn't already declared
    call trace_end ( 'Evaluate_Variable', cond=toggle(gen) )

  end subroutine Evaluate_Variable

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Evaluate_Variable_m

! $Log$
! Revision 2.7  2016/05/25 00:19:30  vsnyder
! Decruftication
!
! Revision 2.6  2014/04/09 00:43:43  vsnyder
! Add Define_Variable_As_String
!
! Revision 2.5  2014/02/27 02:24:44  vsnyder
! Redeclare variable consistently with Declare in tree_checker
!
! Revision 2.4  2014/02/21 19:27:16  vsnyder
! More work on type checking, especially for enumeration-typed variables
!
! Revision 2.3  2014/01/11 01:41:02  vsnyder
! Decruftification
!
! Revision 2.2  2014/01/08 21:09:06  vsnyder
! Add more type checking and tracing
!
! Revision 2.1  2013/10/09 01:06:36  vsnyder
! Initial commit
!
