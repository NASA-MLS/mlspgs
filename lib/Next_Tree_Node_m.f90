! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Next_Tree_Node_m

  implicit NONE
  private

  public :: Init_Next_Tree_Node
  public :: Next_Tree_Node, Next_Tree_Node_State, Traverse_Tree

  type :: Ancestors_t
    integer :: Root = 0    ! of a subtree under consideration
    integer :: Subtree = 0 ! last, not next, subtree of Root that was processed
  end type Ancestors_t

  type :: Next_Tree_Node_State
    type(ancestors_t), allocatable :: Ancestors(:)
  end type Next_Tree_Node_State

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Procedures     =====================================

  subroutine Init_Next_Tree_Node ( State )
    type(next_Tree_Node_State), intent(out) :: State ! Default initialization
  end subroutine Init_Next_Tree_Node

  integer function Next_Tree_Node ( Root, State, Start, NoVariable, TraceLevel )

  ! Get the next tree node, taking IF and CASE constructs into account,
  ! i.e., process them here and select which tree node to return.
  ! To use this:
  !   call Init_Next_Tree_Node ( State )
  !   DO
  !     subtree = next_tree_node ( root, state )
  !     if ( subtree == 0 ) exit ! No more subtrees
  !     process the subtree
  !   END DO

    use Allocate_Deallocate, only: Memory_Units , Test_Allocate, Test_Deallocate
    use Declaration_Table, only: Log_Value, Operator(==), Range, Str_Range, &
      & Value_T
    use Evaluate_Variable_m, only: Evaluate_Variable
    use Expr_M, only: Expr
    use Intrinsic, only: PHYQ_Invalid
    use Toggles, only: Gen, Levels, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    use Tree, only: Node_ID, Nsons, Subtree
    use Tree_Types, only: N_CF, N_Default, N_Else, N_If, N_Select, N_Test, &
      & N_Variable

    integer, intent(in) :: Root    ! Tree index
    type(next_Tree_Node_State), intent(inout) :: State
    integer, intent(in), optional :: Start ! Subtree index of first subtree
    logical, intent(in), optional :: NoVariable ! Don't evaluate variables
    integer, intent(in), optional :: TraceLevel

    logical :: DoVariable          ! Evaluate variables
    integer :: I                   ! Loop inductor
    integer :: Last                ! Last son of Root to process
    integer :: Me = -1             ! String index for trace
    integer :: MyLevel
    integer :: R                   ! Root of subtree of Root
    integer :: Stat                ! From Allocate or Deallocate
    type(ancestors_t), allocatable :: Temp(:)
    integer :: Test_Type, Type
    integer :: Test_Units(2), Units(2)
    double precision :: Test_Value(2), Value(2)
    type(value_t), allocatable :: Test_Values(:), Values(:)

    if ( .not. allocated(state%ancestors) ) then
      allocate ( state%ancestors(1), stat=stat )
      call test_allocate ( stat, moduleName, 'state%ancestors', [1], [1], &
        & storage_size(state%ancestors) )
      state%ancestors(1)%root = root
      state%ancestors(1)%subtree = nsons(root) - 1 ! Assume not N_CF or N_Test
      if ( present(start) ) then
        state%ancestors(1)%subtree = start - 1
      else if ( node_id(root) == n_cf .or. node_id(root) == n_test ) then
        state%ancestors(1)%subtree = 1
      end if
    end if

    doVariable = .true.
    if ( present(noVariable) ) doVariable = .not. noVariable
    myLevel = 0
    if ( present(traceLevel) ) myLevel = traceLevel

    call trace_begin ( me, 'Next_Tree_Node', state%ancestors(1)%root, &
      & index=state%ancestors(1)%subtree+1, &
      & cond=toggle(gen) .and. levels(gen) >= myLevel )
 o: do
      state%ancestors(1)%subtree = state%ancestors(1)%subtree + 1
      last = nsons(state%ancestors(1)%root)
      if ( node_id(state%ancestors(1)%root) == n_cf ) &
        & last = last - 1 ! Skip name at end of section
      if ( state%ancestors(1)%subtree > last ) then
        if ( size(state%ancestors) == 1 ) then
          deallocate ( state%ancestors, stat=stat )
          call test_deallocate ( stat, moduleName, 'state%ancestors', &
            & real(storage_size(state%ancestors)) / memory_units )
          next_Tree_Node = 0
          exit
        end if
        ! Pop the ancestors stack; ascend in the syntax tree.
        allocate ( temp(size(state%ancestors)-1), stat=stat )
        call test_allocate ( stat, moduleName, 'temp', [1], [1], &
          & storage_size(temp) )
        temp = state%ancestors(2:)
        call move_alloc ( temp, state%ancestors )
      else
        next_Tree_Node = subtree(state%ancestors(1)%subtree,state%ancestors(1)%root)
        select case ( node_id(next_tree_node) )
        case ( n_if )
          do i = 1, nsons(next_Tree_Node)
            r = subtree(i,next_Tree_Node)
            select case ( node_id(r) )
            case ( n_test )
              call expr ( subtree(1,r), units, value, type, values=values )
              if ( type /= log_value ) call announce_error ( &
                & subtree(1,r), 'Predicate of IF is not logical type' )
              if ( size(values) /= 1 ) call announce_error ( &
                & subtree(1,r), 'Predicate of IF is not scalar' )
              if ( value(1) /= 0 ) then
                call push
                state%ancestors(1)%subtree=1 ! Start with subtree(2,r)
                cycle o
              end if
            case ( n_else )
              call push
              cycle o
            end select
          end do
        case ( n_select )
          call expr ( subtree(1,next_Tree_Node), units, value, type, values=values )
          if ( size(values) /= 1 ) call announce_error ( &
            & subtree(1,next_Tree_Node), 'Expression in SELECT is not scalar' )
          if ( values(1)%units(1) == range .or. &
             & values(1)%units(1) == str_range) call announce_error ( &
               & subtree(1,next_Tree_Node), 'Expression in SELECT is a range' )
          do i = 2, nsons(next_Tree_Node)
            r = subtree(i,next_Tree_Node)
            select case ( node_id(r) )
            case ( n_test )
              call expr ( subtree(1,r), test_units, test_value, test_type, &
                & values=test_values )
              if ( size(test_values) /= 1 ) call announce_error ( &
                & subtree(1,r), 'Expression in CASE is not scalar' )
              if ( test_values(1)%units(1) == range .or. &
                 & test_values(1)%units(1) == str_range) call announce_error ( &
                   & subtree(1,r), 'Expression in CASE is a range' )
              if ( type /= test_type ) call announce_error ( &
                & subtree(1,r), 'Expressions in SELECT and CASE are not the same type' )
              if ( test_values(1) == values(1) ) then
                call push
                state%ancestors(1)%subtree=1 ! Start with subtree(2,r)
                cycle o
              end if
            case ( n_default )
              call push
              cycle o
            end select
          end do
        case ( n_variable )
          if ( doVariable ) then
            call evaluate_variable ( next_tree_node )
          else
            exit
          end if
        case default ! Not IF or SELECT or VARIABLE
          exit
        end select
      end if
    end do o

    call trace_end ( 'Next_Tree_Node', &
      & cond=toggle(gen) .and. levels(gen) >= myLevel )

  contains

    subroutine Push
      ! Encountered IF with true predicate, ELSE, CASE equal to SELECT,
      ! or CASE DEFAULT; descend into its subtree.
      integer :: N
      n = size(state%ancestors)+1
      allocate ( temp(n), stat=stat )
      call test_allocate ( stat, moduleName, 'temp', [1], [n], &
        & storage_size(temp) )
      temp(2:) = state%ancestors
      call move_alloc ( temp, state%ancestors )
      state%ancestors(1)%root = r
    ! state%ancestors(1)%subtree = 0 ! The default, so we don't need to do it
    end subroutine Push

  end function Next_Tree_Node

  recursive subroutine Traverse_Tree ( Root, Body, Start, NoVariable, &
    & TraceLevel )

    ! Traverse subtrees Start .. nsons(root) of Root.  Handle IF and
    ! SELECT here.  Pass the selected subtree to subroutine "Body."

    use Declaration_Table, only: Log_Value, Operator(==), Range, Str_Range, &
      & Value_T
    use Evaluate_Variable_m, only: Evaluate_Variable
    use Expr_M, only: Expr
    use Intrinsic, only: PHYQ_Invalid
    use Toggles, only: Gen, Levels, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    use Tree, only: Node_ID, Nsons, Subtree
    use Tree_Types, only: N_CF, N_Default, N_Else, N_If, N_Select, N_Test, &
      & N_Variable

    integer, intent(in) :: Root    ! Tree index
    interface
      subroutine Body ( Root )      ! Subroutine to process subtree
        integer, intent(in) :: Root ! Root of subtree to be processed
      end subroutine Body
    end interface
    integer, intent(in), optional :: Start ! Subtree index of first subtree
    logical, intent(in), optional :: NoVariable ! Don't evaluate variables
    integer, intent(in), optional :: TraceLevel

    logical :: DoVariable          ! Evaluate variables
    integer :: First               ! First  son of Root to process
    integer :: I, J                ! Loop inductors
    integer :: Last                ! Last son of Root to process
    integer :: Me = -1             ! String index for trace
    integer :: MyLevel
    integer :: R                   ! Root of subtree of Root
    integer :: Test_Type, Type
    integer :: Test_Units(2), Units(2)
    double precision :: Test_Value(2), Value(2)
    type(value_t), allocatable :: Test_Values(:), Values(:)

    doVariable = .true.
    if ( present(noVariable) ) doVariable = .not. noVariable
    myLevel = 0
    if ( present(traceLevel) ) myLevel = traceLevel

    call trace_begin ( me, 'Traverse_Tree', root, &
      & cond=toggle(gen) .and. levels(gen) >= myLevel )

    if ( present(start) ) then
      first = start
    else
      first = 1 ! Assume not N_CF or N_Test
      if ( node_id(root) == n_cf .or. node_id(root) == n_test ) first = 2
    end if

    last = nsons(root)
    if ( node_id(root) == n_cf ) last = last - 1 ! Skip name at end of section
 o: do i = first, last
      r = subtree(i,root)
      select case ( node_id(r) )
      case ( n_if )
        do j = 1, nsons(r)
          r = subtree(j,r)
          select case ( node_id(r) )
          case ( n_test )
            call expr ( subtree(1,r), units, value, type, values=values )
            if ( type /= log_value ) call announce_error ( &
              & subtree(1,r), 'Predicate of IF is not logical type' )
            if ( size(values) /= 1 ) call announce_error ( &
              & subtree(1,r), 'Predicate of IF is not scalar' )
            if ( value(1) /= 0 ) then
              call body ( subtree(2,r) )
              cycle o
            end if
          case ( n_else )
            call body ( subtree(1,r) )
            cycle o
          end select
        end do
      case ( n_select )
        call expr ( subtree(1,r), units, value, type, values=values )
        if ( size(values) /= 1 ) call announce_error ( &
          & subtree(1,r), 'Expression in SELECT is not scalar' )
        if ( values(1)%units(1) == range .or. &
           & values(1)%units(1) == str_range) call announce_error ( &
             & subtree(1,r), 'Expression in SELECT is a range' )
        do j = 1, nsons(r)
          r = subtree(j,r)
          select case ( node_id(r) )
          case ( n_test )
            call expr ( subtree(1,r), test_units, test_value, test_type, &
              & values=test_values )
            if ( size(test_values) /= 1 ) call announce_error ( &
              & subtree(1,r), 'Expression in CASE is not scalar' )
            if ( test_values(1)%units(1) == range .or. &
               & test_values(1)%units(1) == str_range) call announce_error ( &
                 & subtree(1,r), 'Expression in CASE is a range' )
            if ( type /= test_type ) call announce_error ( &
              & subtree(1,r), 'Expressions in SELECT and CASE are not the same type' )
            if ( test_values(1) == values(1) ) then
              call body ( subtree(2,r) )
              cycle o
            end if
          case ( n_default )
            call body ( subtree(1,r) )
            cycle o
          end select
        end do
      case ( n_variable )
        if ( doVariable ) then
          call evaluate_variable ( r )
        else
          call body ( r )
        end if
      case default ! Not IF or SELECT or VARIABLE
        call body ( r )
      end select
    end do o

    call trace_end ( 'Traverse_Tree', &
      & cond=toggle(gen) .and. levels(gen) >= myLevel )

  end subroutine Traverse_Tree

! ====     Private Procedures     ======================================

  subroutine Announce_Error ( Loc, What )
    use Lexer_Core, only: Print_Source
    use Output_m, only: Output
    use PrintIt_m, only: MLSMSG_Error, PrintItOut
    use Tree, only: Where
    integer, intent(in) :: Loc ! Tree node index
    character(len=*), intent(in) :: What
    call print_source ( where(loc), advance='no', after=': ' )
    ! Print Message, dump calling stack (if any) and stop
    call printitout ( what, MLSMSG_error, exitStatus = 1 )
  end subroutine Announce_Error

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Next_Tree_Node_m

! $Log$
! Revision 2.1  2013/12/12 02:12:07  vsnyder
! Initial commit
!
