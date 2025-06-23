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

  public :: Dump, Init_Next_Tree_Node, IsStateNull
  public :: Next_Tree_Node, Next_Tree_Node_State, Traverse_Tree

  ! Values of Ancestors%what
  integer, parameter :: Go = 0            ! Not immediately within a DO or WHILE
  integer, parameter :: Do_Steps = go + 1 ! Within Do Var := expr, expr [, expr]
  integer, parameter :: Do_Vals = do_Steps + 1 ! Within DO Var := expr
  integer, parameter :: Do_While = do_vals + 1 ! Within DO WHILE ( expr )

  type :: Ancestors_t
    integer :: Root = 0    ! of a subtree under consideration
    integer :: Subtree = 0 ! last subtree of Root that was processed (not the next one)
    integer :: What = go   ! What kind of construct (see parameters above)?
    integer :: Current     ! Current variable value in DO construct
    integer :: Last, Step  ! for DO construct with 2-3 exprs
  end type Ancestors_t

  type :: Next_Tree_Node_State
    type(ancestors_t), allocatable :: Ancestors(:)
  end type Next_Tree_Node_State
  
  interface Dump
    module procedure Dump_ancestor
    module procedure Dump_Next_Tree_Node
  end interface
  
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Procedures     =====================================

  subroutine Dump_ancestor ( ancestor )
    use highOutput, only: startTable, AddRow, AddRow_divider, AddRow_header, &
      & outputTable
    type(Ancestors_t), intent(in) :: ancestor ! What to Dump
    ! Internal variables
    ! Executable
    call startTable
    call addRow_header ( 'ancestor Info', 'c' )
    call addRow_divider ( '-' )
    call addRow ( 'Root     ', ancestor%Root      )
    call addRow ( 'Subtree  ', ancestor%Subtree   )
    call addRow ( 'What     ', ancestor%What      )
    call addRow ( 'Current  ', ancestor%Current   )
    call addRow ( 'Last     ', ancestor%Last      )
    call addRow ( 'Step     ', ancestor%Step      )
    call outputTable ( sep='|', border='-', &
      & Interior=achar(0), headliner='*' )
  end subroutine Dump_ancestor

  subroutine Dump_Next_Tree_Node ( State, Name )
    use highOutput, only: banner
    use Output_m, only: output
    type(next_Tree_Node_State), intent(in) :: State ! What to Dump
    character(len=*), optional, intent(in) :: name
    ! Internal variables
    integer :: i
    character(len=16) :: myName
    ! Executable
    myName = 'state'
    if ( present(name) ) myName = name
    if ( .not. allocated(state%ancestors) ) then
      call output( trim(myname) // ' ancestors not yet allocated', advance='yes' )
      return
    endif
    call banner ( 'Dump of next tree node datatype (' // trim(myName) // ')' )
    do i = 1, size(state%ancestors)
      call Dump_ancestor( state%ancestors(i) )
    enddo
  end subroutine Dump_Next_Tree_Node

  logical function isStateNull ( State )
    type(next_Tree_Node_State), intent(in) :: State ! What to check for being null
    ! Internal variables
    ! Executable
    isStateNull = .not. allocated(state%ancestors)
  end function isStateNull

  subroutine Init_Next_Tree_Node ( State )
    type(next_Tree_Node_State), intent(out) :: State ! Default initialization
  end subroutine Init_Next_Tree_Node

  integer function Next_Tree_Node ( Root, State, Start, NoVariable, TraceLevel )

  ! Get the next tree node, taking IF, SELECT CASE, DO, and WHILE constructs
  ! into account, i.e., process them here and select which tree node to return.
  ! To use this:
  !   call Init_Next_Tree_Node ( State )
  !   DO
  !     subtree = next_tree_node ( root, state )
  !     if ( subtree == 0 ) exit ! No more subtrees
  !     process the subtree
  !   END DO

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    use Declaration_Table, only: Enum_Value, Operator(==), Range, Redeclare, &
      & Str_Range, Value_T, Variable
    use Evaluate_Variable_m, only: Evaluate_Variable
    use Expr_M, only: Expr
    use Intrinsic, only: L_True, T_Boolean  
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSStringLists, only: SwitchDetail
    use Toggles, only: Switches
    use Trace_m, only: Trace_Begin, Trace_End
    use Tree, only: Node_ID, Nsons, Subtree, Sub_Rosa
    use Tree_Types, only: N_CF, N_Cycle, N_Default, N_Do, N_Else, N_Exit, &
      & N_If, N_Named, N_Select, N_Test, N_Variable, N_While

    integer, intent(in) :: Root    ! Tree index
    type(next_Tree_Node_State), intent(inout) :: State
    integer, intent(in), optional :: Start ! Subtree index of first subtree
    logical, intent(in), optional :: NoVariable ! Don't evaluate variables
    integer, intent(in), optional :: TraceLevel

    integer :: Begin               ! Where to begin on a construct, for
                                   ! skipping construct labels
    logical :: DoVariable          ! Evaluate variables
    integer :: I                   ! Loop inductor or subtree index (1 or 2)
    integer :: Last                ! Last son of Root to process
    integer :: Me = -1             ! String index for trace
    integer :: MyLevel             ! Trace level requested by -Snode
    integer :: R                   ! Root of subtree of Root
    integer :: Stat                ! From Allocate or Deallocate
    integer :: SwitchLevel = -2    ! Trace level requested by -Snode
    type(ancestors_t), allocatable :: Temp(:)
    integer :: Test_Type, Type
    integer :: Test_Units(2), Units(2)
    integer :: V                   ! Values subtree of DO expr := values
    double precision :: Test_Value(2), Value(2)
    type(value_t), allocatable :: Test_Values(:), Values(:)

    if ( .not. allocated(state%ancestors) ) then
      allocate ( state%ancestors(1), stat=stat )
      call test_allocate ( stat, moduleName, 'state%ancestors', [1], [1], &
        & storage_size(state%ancestors) / 8 )
      state%ancestors(1)%root = root
      state%ancestors(1)%subtree = nsons(root) - 1 ! Assume not N_CF or N_Test
      if ( present(start) ) then
        state%ancestors(1)%subtree = start - 1
      else if ( node_id(root) == n_cf .or. node_id(root) == n_test ) then
        state%ancestors(1)%subtree = 1
      end if
    end if

    if ( switchLevel < -1 ) switchLevel = switchDetail ( switches, 'node' )
    doVariable = .true.
    if ( present(noVariable) ) doVariable = .not. noVariable
    myLevel = 0
    if ( present(traceLevel) ) myLevel = traceLevel

    call trace_begin ( me, 'Next_Tree_Node', state%ancestors(1)%root, &
      & index=state%ancestors(1)%subtree+1, &
      & cond=switchLevel >= myLevel )
 o: do
      state%ancestors(1)%subtree = state%ancestors(1)%subtree + 1
      r = state%ancestors(1)%root
      last = nsons(r)
      if ( node_id(r) == n_cf ) last = last - 1 ! Skip name at end of section
      if ( state%ancestors(1)%subtree > last ) then
        select case ( state%ancestors(1)%what )
        case ( go ) ! Not a looping construct, just finish it
          if ( pop() ) exit o
        case ( do_steps ) ! Check whether we've done all the steps
          state%ancestors(1)%current = state%ancestors(1)%current + &
                                       state%ancestors(1)%step
          if ( ( state%ancestors(1)%current - state%ancestors(1)%last ) * &
               & sign(1,state%ancestors(1)%step) > 0 ) then
          if ( pop() ) exit o
          else
            i = 1
            if ( node_id(subtree(1,state%ancestors(1)%root)) == n_named ) i = 2
            call set_variable_numeric ( subtree(1,subtree(i,state%ancestors(1)%root)), &
              & dble(state%ancestors(1)%current) )
            state%ancestors(1)%subtree = i ! Do the block again
          end if
        case ( do_vals ) ! Check whether we've used all the values
          i = 1
          if ( node_id(subtree(1,state%ancestors(1)%root)) == n_named ) i = 2
          ! node_id(v) == n_array, because it has to be something, but
          ! who cares?
          v = subtree(2,subtree(i,r))
          ! %Current is the index of the last expression evaluated
          state%ancestors(1)%current = state%ancestors(1)%current + 1
          if ( state%ancestors(1)%current > nsons(v) ) then ! Used them all
          if ( pop() ) exit o
          else ! Use the next one
            call expr ( subtree(state%ancestors(1)%current,v), units, value, &
              & type, values=values )
            call redeclare ( sub_rosa(subtree(1,subtree(i,r))), value(1), &
              & variable, type, values=values )
            state%ancestors(1)%subtree = i ! Do the block again
          end if
        case ( do_while ) ! Check whether the expr is true
          i = 1
          if ( node_id(subtree(1,r)) == n_named ) i = 2
          ! Type checker has verified that expression type is boolean
          call expr ( subtree(i,r), units, value )
          if ( nint(value(1)) /= l_true ) then ! Value is false, finish the construct
          if ( pop() ) exit o
          else
            state%ancestors(1)%subtree = i ! Do the block again
          end if
        end select
      else
        next_Tree_Node = subtree(state%ancestors(1)%subtree,state%ancestors(1)%root)
        select case ( node_id(next_tree_node) )
        case ( n_cycle )
          if ( cycle_or_exit () ) exit o
        case ( n_do )
          call push ( next_tree_node )
          begin = 1
          if ( node_id(subtree(1,next_tree_node)) == n_named ) begin = 2
          r = subtree(begin,next_tree_node) ! N_Do_Head vertex
          if ( nsons(r) == 2 ) then
            ! Variable := value (maybe an array)
            state%ancestors(1)%what = do_vals
            ! Pretend we've done its block once, so we check the count at the top
            state%ancestors(1)%subtree = nsons(next_tree_node)
            state%ancestors(1)%current = 0 ! Which value we just finished
          else
            ! Variable := expr, expr [, expr]
            state%ancestors(1)%what = do_steps
            ! Type checker has verified that exprs are numeric
            call expr ( subtree(2,r), units, value, values=values )
            state%ancestors(1)%current = nint(value(1))
            call set_variable_numeric ( subtree(1,r), value(1) )
            call expr ( subtree(3,r), units, value, values=values )
            state%ancestors(1)%last = nint(value(1))
            value(1) = 1
            if ( nsons(r) > 3 ) &
              & call expr ( subtree(4,r), units, value, values=values )
            state%ancestors(1)%step = nint(value(1))
            if ( state%ancestors(1)%step == 0 ) call MLSMessage ( MLSMSG_Error, &
              & moduleName, 'STEP in a DO construct is zero' )
            ! Pretend we've done its block once, so we check the count at the top
            state%ancestors(1)%subtree = nsons(next_tree_node)
            ! Current will be bumped at the top
            state%ancestors(1)%current = state%ancestors(1)%current - &
              state%ancestors(1)%step
          end if
        case ( n_exit )
          if ( cycle_or_exit () ) exit o
          if ( pop() ) exit o
        case ( n_if )
          begin = 1
          if ( node_id(subtree(1,next_tree_node)) == n_named ) begin = 2
          do i = begin, nsons(next_Tree_Node)
            r = subtree(i,next_Tree_Node)
            select case ( node_id(r) )
            case ( n_test )
              ! Type checked in tree_checker
              call expr ( subtree(1,r), units, value, values=values )
              if ( values(1)%what /= enum_value .or. &
                 & values(1)%type /= t_boolean ) call announce_error ( &
                & subtree(1,r), 'Predicate of IF is not logical' )
              if ( size(values) /= 1 ) call announce_error ( &
                & subtree(1,r), 'Predicate of IF is not scalar' )
              if ( nint(value(1)) == l_true ) then
                call push ( r )
                state%ancestors(1)%subtree=1 ! Start with subtree(2,r)
                cycle o
              end if
            case ( n_else )
              call push ( r )
              cycle o
            end select
          end do
        case ( n_select )
          begin = 1
          if ( node_id(subtree(1,next_tree_node)) == n_named ) begin = 2
          call expr ( subtree(begin,next_Tree_Node), units, value, type, values=values )
          if ( size(values) /= 1 ) call announce_error ( &
            & subtree(begin,next_Tree_Node), 'Expression in SELECT is not scalar' )
          if ( values(1)%units(1) == range .or. &
             & values(1)%units(1) == str_range) call announce_error ( &
               & subtree(begin,next_Tree_Node), 'Expression in SELECT is a range' )
          do i = begin+1, nsons(next_Tree_Node)
            r = subtree(i,next_Tree_Node)
            select case ( node_id(r) )
            case ( n_test )
              call expr ( subtree(1,r), test_units, test_value, test_type, &
                & values=test_values )
              if ( size(test_values) /= 1 ) call announce_error ( &
                & subtree(1,r), 'Expression in CASE is not scalar' )
              if ( test_type /= type ) call announce_error ( subtree(1,r), &
                & 'Type of expression in CASE is not the same as in SELECT' )
              if ( test_values(1)%units(1) == range .or. &
                 & test_values(1)%units(1) == str_range) call announce_error ( &
                   & subtree(1,r), 'Expression in CASE is a range' )
              if ( test_values(1) == values(1) ) then
                call push ( r )
                state%ancestors(1)%subtree=1 ! Start with subtree(2,r)
                cycle o
              end if
            case ( n_default )
              call push ( r )
              cycle o
            end select
          end do
        case ( n_variable )
          if ( doVariable ) then
            call evaluate_variable ( next_tree_node )
          else
            exit
          end if
        case ( n_while )
          call push ( next_tree_node )
          state%ancestors(1)%what = do_while
          ! Pretend we've done its block once, so we check the expr at the top
          state%ancestors(1)%subtree = nsons(next_tree_node)
        case default ! Not DO, IF, SELECT, VARIABLE, or WHILE
          exit
        end select
      end if
    end do o

    call trace_end ( 'Next_Tree_Node', index=next_tree_node, &
      & cond=switchLevel >= myLevel )

  contains

    logical function Cycle_Or_Exit ()
      ! Find the referenced construct.  If the ancestors stack goes empty,
      ! the result variable is true.
      cycle_or_exit = .false.
      do
        if ( state%ancestors(1)%what /= go ) then
          r = state%ancestors(1)%root
          state%ancestors(1)%subtree = nsons(r) ! Done with the block
          if ( nsons(next_tree_node) == 0 ) exit ! bare cycle
          if ( node_id(subtree(1,r)) == n_named ) then
            if ( sub_rosa(subtree(1,next_tree_node)) == &
               & sub_rosa(subtree(1,subtree(1,r))) ) exit ! labels match
          end if
        end if
        cycle_or_exit = pop()
        if ( cycle_or_exit ) exit
      end do
    end function Cycle_Or_Exit

    logical function Pop ()
      ! Pop the ancestors stack.  The result variable is true, and
      ! next_tree_node = 0, if it goes empty.
      pop = .false.
      if ( size(state%ancestors) == 1 ) then
        deallocate ( state%ancestors, stat=stat )
        call test_deallocate ( stat, moduleName, 'state%ancestors', &
          & storage_size(state%ancestors) / 8 )
        next_tree_node = 0
        pop = .true.
      else
        ! Pop the ancestors stack; ascend in the syntax tree.
        allocate ( temp(size(state%ancestors)-1), stat=stat )
        call test_allocate ( stat, moduleName, 'temp', [1], [1], &
          & storage_size(temp) / 8 )
        temp = state%ancestors(2:)
        call move_alloc ( temp, state%ancestors )
      end if
    end function Pop

    subroutine Push ( Root )
      ! Encountered IF with true predicate, ELSE, CASE equal to SELECT,
      ! or CASE DEFAULT; descend into its subtree.
      integer, intent(in) :: Root
      integer :: N
      n = size(state%ancestors)+1
      allocate ( temp(n), stat=stat )
      call test_allocate ( stat, moduleName, 'temp', [1], [n], &
        & storage_size(temp) / 8 )
      temp(2:) = state%ancestors
      call move_alloc ( temp, state%ancestors )
      state%ancestors(1)%root = root
    ! state%ancestors(1)%subtree = 0 ! The default, so we don't need to do it
    ! state%ancestors(1)%what = go   ! The default, so we don't need to do it
    end subroutine Push

  end function Next_Tree_Node

  subroutine Set_Variable_Numeric ( Root, Value )

    use Declaration_Table, only: Num_Value, Redeclare, Value_Allocate, &
      & Value_T, Variable
    use Intrinsic, only: PHYQ_Dimensionless
    use Tree, only: Sub_Rosa

    integer, intent(in) :: Root ! Name node of the variable
    double precision, intent(in) :: Value
    type(value_t), allocatable :: Values(:)
    call value_allocate ( values, 1, 'Values', moduleName )
    values(1) = value_t(num_value,num_value,[value,0.0d0], &
                       &[PHYQ_Dimensionless,PHYQ_Dimensionless],0)
    call redeclare ( sub_rosa(root), value, variable, num_value, root, &
      & values=values )

  end subroutine Set_Variable_Numeric

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                                !!!!!
!!!!!  The recursive iterator is obsolete!  It doesn't handle DO     !!!!!
!!!!!  constructs, WHILE constructs, CYCLE statements, or EXIT       !!!!!
!!!!!  statements!  It handles boolean type the old way (log_value)! !!!!!
!!!!!                                                                !!!!!
!!!!!  Fortunately, it's not used.  Its purpose is to outline how    !!!!!
!!!!!  it might be done if we ever want to pass the subtree          !!!!!
!!!!!  handler as a dummy procedure.                                 !!!!!
!!!!!                                                                !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  recursive subroutine Traverse_Tree ( Root, Body, Start, NoVariable, &
    & TraceLevel )

    ! Traverse subtrees Start .. nsons(root) of Root.  Handle IF and
    ! SELECT here.  Pass the selected subtree to subroutine "Body."

    use Declaration_Table, only: Log_Value, Operator(==), Range, Str_Range, &
      & Value_T
    use Evaluate_Variable_m, only: Evaluate_Variable
    use Expr_M, only: Expr
    use MLSStringLists, only: SwitchDetail
    use Toggles, only: Switches
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
    integer :: SwitchLevel = -2    ! Trace level requested by -Snode
    integer :: Test_Type, Type
    integer :: Test_Units(2), Units(2)
    double precision :: Test_Value(2), Value(2)
    type(value_t), allocatable :: Test_Values(:), Values(:)

    if ( switchLevel < -1 ) switchLevel = switchDetail ( switches, 'node' )
    doVariable = .true.
    if ( present(noVariable) ) doVariable = .not. noVariable
    myLevel = 0
    if ( present(traceLevel) ) myLevel = traceLevel

    call trace_begin ( me, 'Traverse_Tree', root, &
      & cond=switchLevel >= myLevel )

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
      & cond=switchLevel >= myLevel )

  end subroutine Traverse_Tree

! ====     Private Procedures     ======================================

  subroutine Announce_Error ( Loc, What )
    use Lexer_Core, only: Print_Source
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
! Revision 2.9  2017/01/19 23:34:31  pwagner
! Improve appeearance when dumping ancestor
!
! Revision 2.8  2016/05/25 00:19:30  vsnyder
! Decruftication
!
! Revision 2.7  2015/07/31 20:41:29  pwagner
! Added Dump, isStateNull
!
! Revision 2.6  2014/09/05 18:37:00  vsnyder
! Measure memory usage in bytes instead of Memory_Units
!
! Revision 2.5  2014/02/28 00:02:11  vsnyder
! Turn off type checking already done by type checker.  Check that expr
! in CASE is the same type as in SELECT (which the type checker can't do).
!
! Revision 2.4  2014/02/27 02:25:57  vsnyder
! Handle labels on IF, and SELECT CASE constructs
!
! Revision 2.3  2014/02/21 19:26:30  vsnyder
! Add CYCLE, DO, EXIT, WHILE
!
! Revision 2.2  2014/01/11 01:41:02  vsnyder
! Decruftification
!
! Revision 2.1  2013/12/12 02:12:07  vsnyder
! Initial commit
!
