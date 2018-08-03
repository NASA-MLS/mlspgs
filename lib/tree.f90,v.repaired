! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module TREE

  use error_handler, only: compiler, error_intro
  use machine, only: io_error
  use printIt_m, only: printItOut, MLSMSG_Error
  use output_m, only: newline, output
  use string_table, only: display_string, get_string, lookup_and_insert
  use symbol_table, only: set_symbol, symbol
  use symbol_types, only: treeNode
  use toggles, only: con, toggle
  use tree_types, only: first_tree_node, last_tree_node, n_null, tree_init, &
    & tree_map
  implicit NONE
  private

  public :: ADD_SONS_FROM_STACK, ALLOCATE_TREE, BUILD_TREE, COPY_TO_STACK
  public :: DEALLOCATE_TREE, DECORATE, DECORATION, DECORATION_TX_TX, DELETE_TREE
  public :: DELETE_TREE_STACK, DUMP_STACK, DUMP_TOP_STACK, DUMP_TOP_STACK_NAME
  public :: DUMP_TREE_NODE, DUMP_TREE_NODE_NAME, GET_TREE_NODE_NAME, INIT_TREE
  public :: INSERT_NODE, NODE_ID, NODE_IN_TREE, NODE_KIND, NSONS, POP, PRINT_SUBTREE
  public :: PUSH_PSEUDO_TERMINAL, REPLACE_SONS, SOURCE_REF, STACK_FILE
  public :: STACK_SOURCE_REF,STACK_SUB_ROSA, STACK_SUBTREE, STACK_SUBTREE_TX
  public :: SUB_ROSA, SUBTREE, THE_FILE, TREE_NODE_NAME,TREE_TEXT, TX, WHERE

  ! Tree node kinds:
  integer, public, parameter :: PSEUDO = 0   ! Tree node is pseudo terminal
  integer, public, parameter :: INTERNAL = 1 ! Tree node has sons
  integer, public, parameter :: MORE = 2     ! Used to link tree after
                                             !  transformations
  integer, public, parameter :: EMPTY = 3    ! Empty node

  ! How many nodes in the tree stack?
  integer, public, save :: N_TREE_STACK
  ! Index of the "null tree node"
  integer, public, parameter :: NULL_TREE = 0

  interface Add_Sons_From_Stack
    module procedure Add_Sons_From_Stack_I, Add_Sons_From_Stack_TX
  end interface

  interface Copy_To_Stack
    module procedure Copy_To_Stack_I, Copy_To_Stack_TX
  end interface

  interface Decorate
    module procedure Decorate_I, Decorate_TX, Decorate_TX_TX
  end interface

  interface Decoration
    module procedure Decoration_I, Decoration_TX
  end interface

  interface Delete_Tree
    module procedure Delete_Tree_I, Delete_Tree_TX
  end interface

  interface Dump_Tree_Node
    module procedure Dump_Tree_Node_I, Dump_Tree_Node_TX
  end interface

  interface Dump_Tree_Node_Name
    module procedure Dump_Tree_Node_Name_I, Dump_Tree_Node_Name_TX
  end interface

  interface Get_Tree_Node_Name
    module procedure Get_Tree_Node_Name_I, Get_Tree_Node_Name_TX
  end interface

  interface Insert_Node
    module procedure Insert_Node_I, Insert_Node_TX
  end interface

  interface Node_ID
    module procedure Node_ID_I, Node_ID_TX
  end interface

  interface Node_Kind
    module procedure Node_Kind_I, Node_Kind_TX
  end interface

  interface Nsons
    module procedure Nsons_I, Nsons_TX
  end interface

  interface Print_Subtree
    module procedure Print_Subtree_I, Print_Subtree_TX
  end interface

  interface Push_Pseudo_Terminal
    module procedure Push_Pseudo_Terminal_Integer, Push_Pseudo_Terminal_Where
  end interface

  interface Replace_Sons
    module procedure Replace_Sons_I, Replace_Sons_TX
  end interface

  interface Source_Ref
    module procedure Source_Ref_I, Source_Ref_TX
  end interface

  interface Sub_Rosa
    module procedure Sub_Rosa_I, Sub_Rosa_TX
  end interface

  interface Subtree
    module procedure Subtree_I, Subtree_TX
  end interface

  interface The_File
    module procedure The_File_I, The_File_TX
  end interface

  interface Tree_Node_Name
    module procedure Tree_Node_Name_I, Tree_Node_Name_TX
  end interface

  interface Where
    module procedure Where_I, Where_TX
  end interface

  type :: TX              ! For tree node index, to give strong typing
    integer :: I = 0      ! The real tree node index
  end type TX

  type :: TREE_NODE
    integer :: NODE       ! What kind of tree node
    integer :: DECOR      ! Decoration, an integer or tree node index
    integer :: FILE = 0   ! String index of file
    integer :: NSONS      ! How many sons
    integer :: SOURCE = 0 ! 256*line + column
    integer :: KIND       ! PSEUDO, INTERNAL, MORE or OTHER
    integer :: LINK       ! Sub_rosa if PSEUDO, Left_son if INTERNAL, next
                          !  if MORE, not used if EMPTY
  end type TREE_NODE

  type(TREE_NODE), save, allocatable :: THE_TREE(:)
  integer, save :: TREE_POINT ! Position in tree of last node added.
  integer, save :: TREE_SP    ! Next available space in the orchard stack.

  ! Parameters for tree error codes.  See TREE_ERROR
  integer, parameter :: IS_PSEUDO = 1
  integer, parameter :: NO_MORE_SONS = IS_PSEUDO + 1
  integer, parameter :: NOT_PSEUDO = NO_MORE_SONS + 1
  integer, parameter :: NO_TREE_SPACE = NOT_PSEUDO + 1
  integer, parameter :: UNDERFLOW = NO_TREE_SPACE + 1

  integer, save, private :: TREE_TEXTS(FIRST_TREE_NODE:LAST_TREE_NODE)

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine ADD_SONS_FROM_STACK_I ( T, N )
  ! Add the top N tree nodes on the stack as new sons of T, after the last
  ! son already present.  Pop the N nodes from the stack.
    integer, intent(in) :: T            ! The tree node to get the new sons
    integer, intent(in) :: N            ! How many sons to get
    integer :: I                        ! Loop inductor
    logical :: IN_PLACE                 ! "new nodes can be stored adjacent
                                        !  to old ones"
    integer :: NS                       ! NSONS(T)
    integer :: S                        ! Temp
    if ( tree_point >= tree_sp ) then
      call tree_error ( no_tree_space, null_tree )
!     call double_tree
    end if
    ns = nsons(t)
    s = subtree(ns,t)
    if ( s == tree_point ) then
      call pop_to_tree ( n )
      call set_nsons ( t, ns+n )
    else
    ! determine whether the next "n" cells are "empty".  We assume
    ! we're not going to run off the end of the tree -- Replace_Sons
    ! should decrease tree_point if it would make empties at the end of
    ! the tree space.
      in_place = .true.
      i = 1
      do
        if ( i > n ) then
      exit
        end if
        s = s + 1
        i = i + 1
        if ( the_tree(s)%kind /= empty ) then
          in_place = .false.
      exit
        end if
      end do
      s = subtree(ns,t)
      if ( in_place ) then
      ! there are "n" "empty" cells adjacent to the tree -- use them
        do i = 1, n
          s = s + 1
          the_tree(s) = the_tree(tree_sp+i-1)
        end do
        n_tree_stack = n_tree_stack - n
        tree_sp = tree_sp + n
        call set_nsons ( t, ns+n )
      else
      ! move the last son to the end of the tree space
        tree_point = tree_point + 1
        the_tree(tree_point) = the_tree(s)
        ! change the last son to a "more" node
        the_tree(s)%kind = more
        the_tree(s)%link  = tree_point
        the_tree(s)%nsons  = n+1
        ! pop the additional sons from the stack
        call pop_to_tree ( n )
      end if ! in_place
    end if
  end subroutine ADD_SONS_FROM_STACK_I

  subroutine ADD_SONS_FROM_STACK_TX ( T, N )
  ! Add the top N tree nodes on the stack as new sons of T, after the last
  ! son already present.  Pop the N nodes from the stack.
    type(tx), intent(in) :: T           ! The tree node to get the new sons
    integer, intent(in) :: N            ! How many sons to get
    call add_sons_from_stack ( t%i, n )
  end subroutine ADD_SONS_FROM_STACK_TX

  subroutine ALLOCATE_TREE ( N_TREE, STATUS )
  ! Allocate a TREE array with N_TREE elements
    integer, intent(in) :: N_TREE
    integer, intent(out), optional :: STATUS ! from ALLOCATE
    integer :: STAT
    if ( allocated(the_tree) ) deallocate ( the_tree )
    allocate ( the_tree(0:n_tree), stat=stat )
    if ( stat /= 0 ) then
      if ( present(status) ) then
        status = stat
        return
      end if
      call io_error ( 'TREE%ALLOCATE_TREE-E- Unable to allocate tree', &
        stat, ' ' )
      call PrintItOut ( 'TREE%ALLOCATE_TREE-E- Unable to allocate tree', &
        & 1, exitStatus=1 )
    end if
    call init_tree
  end subroutine ALLOCATE_TREE

  subroutine BUILD_TREE ( NEW_NODE, NSONS, DECORATION, Trace )
  ! Pop NSONS nodes from the orchard stack to the tree.  Build a new tree
  ! node of type NEW_NODE over them.  Leave the new node on the stack at
  ! TREE_SP+1.
    integer, intent(in) :: NEW_NODE
    integer, intent(in) :: NSONS
    integer, intent(in), optional :: DECORATION
    integer, intent(in), optional :: Trace
    integer :: FILE, LEFT_SON, MY_DECOR, SOURCE_REF
    my_decor = null_tree
    if ( present(decoration) ) my_decor = decoration
    if ( tree_sp + nsons > ubound(the_tree,1) ) then
      call tree_error ( underflow, null_tree )
    end if
    if ( nsons == 0 ) then; file = 0; source_ref = 0
    else
      file = the_tree(tree_sp+nsons)%file
      source_ref = the_tree(tree_sp+nsons)%source
    end if
    left_son = tree_point + 1
    call pop_to_tree ( nsons )
    ! Push new node onto the stack
    if ( tree_sp <= tree_point ) then
      call tree_error ( no_tree_space, null_tree )
!     call double_tree
    end if
                                   ! node      decor   file, nsons  source
    the_tree(tree_sp) = tree_node( new_node, my_decor, file, nsons, source_ref, &
                                 !    kind   left_son
                                   internal, left_son )
    if ( present(trace) ) then
      if ( trace > 0 ) then
        call output ( 'Build ' )
        call print_subtree ( tree_sp, 0, .true. )
      end if
    end if
    tree_sp = tree_sp - 1     ! Stack push
    n_tree_stack = n_tree_stack + 1 - nsons
  end subroutine BUILD_TREE

  subroutine COPY_TO_STACK_I ( T, K, M )
  ! Copy sons k..m of t to the top of the stack
    integer, intent(in) :: T            ! Parent of nodes to be copied
    integer, intent(in) :: K            ! First son to be copied
    integer, intent(in) :: M            ! Last son to be copied
    integer :: I                        ! Loop inductor
    integer :: S                        ! SUBTREE(I,T)
    do i = k, m
      s = subtree( i, t )
      the_tree(tree_sp) = the_tree(s)
      tree_sp = tree_sp - 1             ! Push stack
      n_tree_stack = n_tree_stack + 1
    end do
  end  subroutine COPY_TO_STACK_I

  subroutine COPY_TO_STACK_TX ( T, K, M )
  ! Copy sons k..m of t to the top of the stack
    type(tx), intent(in) :: T           ! Parent of nodes to be copied
    integer, intent(in) :: K            ! First son to be copied
    integer, intent(in) :: M            ! Last son to be copied
    call copy_to_stack ( t%i, k, m )
  end subroutine COPY_TO_STACK_TX

  subroutine DEALLOCATE_TREE
    deallocate( the_tree )
  end subroutine DEALLOCATE_TREE

  subroutine DECORATE_I ( WHERE, DECORATION )
  ! Decorate tree node at WHERE with DECORATION
    integer, intent(in) :: WHERE
    integer, intent(in) :: DECORATION
    the_tree(where) % decor = decoration
    if ( toggle(con) ) then
    ! in Fortran 2000:
    ! write ( prunit, '("Decorate ", i0, ": ")', advance="no" ) where
      call output ( "Decorate " )
      call output ( where )
      call output ( ': ' )
      call dump_tree_node ( where, 0 )
    ! write ( prunit, '(" (", i0, ")")' ) decoration ! in Fortran 2000
      call output ( ' (' )
      call output ( decoration )
      call output ( ')', advance='yes' )
    end if
  end subroutine DECORATE_I

  subroutine DECORATE_TX ( WHERE, DECORATION )
  ! Decorate tree node at WHERE with DECORATION
    type(tx), intent(in) :: WHERE
    integer, intent(in) :: DECORATION
    call decorate ( where%i, decoration )
  end subroutine DECORATE_TX

  subroutine DECORATE_TX_TX ( WHERE, DECORATION )
  ! Decorate tree node at WHERE with DECORATION
    type(tx), intent(in) :: WHERE
    type(tx), intent(in) :: DECORATION
    call decorate ( where%i, decoration%i )
  end subroutine DECORATE_TX_TX

  pure integer function DECORATION_I ( WHERE ) result ( Decoration )
  ! Return the decoration of the tree node at WHERE
    integer, intent(in) :: WHERE
    decoration = the_tree(where) % decor
  end function DECORATION_I

  pure integer function DECORATION_TX ( WHERE ) result ( Decoration )
  ! Return the decoration of the tree node at WHERE
    type(tx), intent(in) :: WHERE
    decoration = the_tree(where%i) % decor
  end function DECORATION_TX

  pure type(tx) function DECORATION_TX_TX ( WHERE ) result ( Decoration )
  ! Return the decoration of the tree node at WHERE
    type(tx), intent(in) :: WHERE
    decoration%i = the_tree(where%i) % decor
  end function DECORATION_TX_TX

  subroutine DELETE_TREE_I ( WHERE )
  ! Discard everything at WHERE and above in the tree
    integer, intent(in) :: WHERE
    tree_point = where - 1
  end subroutine DELETE_TREE_I

  subroutine DELETE_TREE_TX ( WHERE )
  ! Discard everything at WHERE and above in the tree
    type(tx), intent(in) :: WHERE
    tree_point = where%i - 1
  end subroutine DELETE_TREE_TX

  subroutine DELETE_TREE_STACK
  ! Discard everything in the tree stack
    tree_sp = ubound(the_tree,1)
    n_tree_stack = 0
  end subroutine DELETE_TREE_STACK

  subroutine DUMP_STACK ( Subtrees )
    ! Dump the tree stack.  If Subtrees is present and true, dump trees
    ! rooted in the stack.
    logical, intent(in), optional :: Subtrees
    integer :: Depth, I
    logical :: My_Trees
    my_trees = .false.
    if ( present(subtrees) ) my_trees = subtrees
    do i = tree_sp+1, ubound(the_tree,1)
      if ( my_trees ) then
        depth = 0
        call print_subtree ( i, depth )
      else
        call output ( i, 5 )
        call output ( ':' )
        call dump_tree_node ( i, 0, advance='yes' )
      end if
    end do
  end subroutine DUMP_STACK

  subroutine DUMP_TOP_STACK ( INDENT, ADVANCE )
    ! Indent INDENT spaces, then dump the top node of the tree stack
    integer, intent(in) :: INDENT
    character(len=*), intent(in), optional :: ADVANCE
    call dump_tree_node ( tree_sp+1, indent, advance )
  end subroutine DUMP_TOP_STACK

  subroutine DUMP_TOP_STACK_NAME ( ADVANCE, BEFORE )
    ! Dump the name of the top node of the tree stack
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: BEFORE
    call display_string ( tree_texts(the_tree(tree_sp+1) % node), &
                        & advance=advance, before=before )
  end subroutine DUMP_TOP_STACK_NAME

  subroutine DUMP_TREE_NODE_I ( WHERE, INDENT, ADVANCE, TYPE_NAME )
  ! Indent INDENT spaces, then dump the tree node at WHERE
    integer, intent(in) :: WHERE
    integer, intent(in) :: INDENT
    character(len=*), intent(in), optional :: ADVANCE
    optional :: Type_Name
    interface
      integer function Type_Name ( Decor )
      ! Return the string index to print for the decoration
        integer, intent(in) :: Decor ! Tree(where)%Decor
      end function Type_Name
    end interface
    integer :: I
    do i = 1, indent; call output ( '.' ); end do
    select case ( the_tree(where) % kind )
    case ( empty )
      call output ( 'empty' )
    case ( more )
      call output ( the_tree(where) % nsons, before='more, nsons = ' )
      call output ( the_tree(where) % link, before=', next = ' )
    case ( pseudo, internal )
      call dump_tree_node_name ( where )
      if ( the_tree(where) % kind == pseudo ) then
        call display_string ( sub_rosa(where), before=' ' )
      else
        call output ( the_tree(where) % nsons, before=', ' )
        call output ( ' sons' )
      end if
      if ( the_tree(where)%decor /= null_tree ) then
        call output ( the_tree(where) % decor, before=' decor=' )
        i = 0
        if ( present(type_name) .and. the_tree(where)%kind == internal ) &
          & i = type_name(the_tree(where)%decor)
        if ( i /= 0 ) call display_string ( i, before=' type=' )
      end if
    end select
    if ( the_tree(where)%source /= 0 ) then
      call output ( the_tree(where)%source/256, before=' line ' )
      call output ( mod(the_tree(where)%source,256), before=' column ' )
    end if
    if ( the_tree(where)%file /= 0 ) &
      & call display_string ( the_tree(where)%file, before=' in ' )
    call output ( '', advance=advance )
  end subroutine DUMP_TREE_NODE_I

  subroutine DUMP_TREE_NODE_TX ( WHERE, INDENT, ADVANCE, TYPE_NAME )
  ! Indent INDENT spaces, then dump the tree node at WHERE
    type(tx), intent(in) :: WHERE
    integer, intent(in) :: INDENT
    character(len=*), intent(in), optional :: ADVANCE
    optional :: TYPE_NAME
    interface
      integer function Type_Name ( Decor )
      ! Return the string index to print for the decoration
        integer, intent(in) :: Decor ! Tree(where)%Decor
      end function Type_Name
    end interface
    call dump_tree_node ( where%i, indent, advance, type_name=type_name )
  end subroutine DUMP_TREE_NODE_TX

  subroutine DUMP_TREE_NODE_NAME_I ( WHERE, ADVANCE, BEFORE )
  ! Dump the name of the tree node at WHERE
    integer, intent(in) :: WHERE
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: BEFORE
    call display_string ( tree_texts(the_tree(where) % node), advance=advance, &
                        & before=before )
  end subroutine DUMP_TREE_NODE_NAME_I

  subroutine DUMP_TREE_NODE_NAME_TX ( WHERE, ADVANCE, BEFORE )
  ! Dump the name of the tree node at WHERE
    type(tx), intent(in) :: WHERE
    character(len=*), intent(in), optional :: ADVANCE
    character(len=*), intent(in), optional :: BEFORE
    call display_string ( tree_texts(the_tree(where%i) % node), advance=advance, &
                        & before=before )
  end subroutine DUMP_TREE_NODE_NAME_TX

  subroutine GET_TREE_NODE_NAME_I ( WHERE, STRING )
  ! Get the name of the tree node at WHERE
    integer, intent(in) :: WHERE
    character(len=*), intent(out), optional :: STRING
    call get_string ( tree_texts(the_tree(where) % node), string )
  end subroutine GET_TREE_NODE_NAME_I

  subroutine GET_TREE_NODE_NAME_TX ( WHERE, STRING )
  ! Get the name of the tree node at WHERE
    type(tx), intent(in) :: WHERE
    character(len=*), intent(out), optional :: STRING
    call get_string ( tree_texts(the_tree(where%i) % node), string )
  end subroutine GET_TREE_NODE_NAME_TX

  subroutine INIT_TREE
    logical :: FOUND     ! Did lookup_and_insert find it?
    integer :: I         ! Loop inductor
    integer :: WHERE     ! Where did lookup_and_insert find it?
    do i = first_tree_node, last_tree_node
      call tree_init (i)
      call lookup_and_insert ( where, found, .false. )
      ! It's OK if it found one -- maybe it's a terminal text, too.
      if ( .not. found ) call set_symbol ( where, treeNode )
      tree_texts(i) = where
    end do
    tree_point = null_tree
    call delete_tree_stack
                                      ! node  decor   file nsons source kind
    the_tree(tree_point) = tree_node( n_null, null_tree, 0,   0, 0, internal, &
                                      ! left_son
                                      null_tree )
  end subroutine INIT_TREE

  subroutine INSERT_NODE_I ( NEW_NODE, T, K, M )
  ! Insert a node having ID = newNode as the k'th son of t.  Make the sons
  ! k..m of t sons of newNode.  Reduce the number of sons of t to
  ! n-k+m-1.  Don't use this procedure during parsing:  It checks nTree,
  ! not treesp!
    integer, intent(in) :: NEW_NODE     ! ID of the node to be inserted
    integer, intent(in) :: T            ! The parent of the inserted node
    integer, intent(in) :: K            ! Which son of T NEW_NODE is to be
    integer, intent(in) :: M            ! Sons k..m of T become sons of
                                        !  NEW_NODE
    integer :: N                        ! NSONS(T)
    n = nsons(t)
    if ( k <= 0 .or. k >= m .or. m >= n ) then
      call tree_error ( no_more_sons, t )
    end if
    call copy_to_stack ( t, k, m )
    n = m - k + 1
    call build_tree ( new_node, n )
    call replace_sons ( t, m, k, tree_sp + 1 )
    call pop ( 1 )
  end subroutine INSERT_NODE_I

  subroutine INSERT_NODE_TX ( NEW_NODE, T, K, M )
  ! Insert a node having ID = newNode as the k'th son of t.  Make the sons
  ! k..m of t sons of newNode.  Reduce the number of sons of t to
  ! n-k+m-1.  Don't use this procedure during parsing:  It checks nTree,
  ! not treesp!
    integer, intent(in) :: NEW_NODE     ! ID of the node to be inserted
    type(tx), intent(in) :: T           ! The parent of the inserted node
    integer, intent(in) :: K            ! Which son of T NEW_NODE is to be
    integer, intent(in) :: M            ! Sons k..m of T become sons of
                                        !  NEW_NODE
    call insert_node ( new_node, t%i, k, m )
  end subroutine INSERT_NODE_TX

  pure integer function NODE_ID_I ( WHERE ) result ( Node_ID )
  ! Return the node id of the tree node at WHERE
    integer, intent(in) :: WHERE
    node_id = the_tree(where) % node
  end function NODE_ID_I

  pure integer function NODE_ID_TX ( WHERE ) result ( Node_ID )
  ! Return the node id of the tree node at WHERE
    type(tx), intent(in) :: WHERE
    node_id = the_tree(where%i) % node
  end function NODE_ID_TX

  function Node_in_tree ( Tree ) result ( inside )
  ! Return TRUE if the node num is a valid index into the tree array
    integer, intent(in) :: Tree ! Tree node index tested
    logical :: inside
    inside = ( Tree >= lbound(the_tree, 1) .and. Tree <= ubound(the_tree, 1) )
  end function Node_in_tree

  pure integer function NODE_KIND_I ( WHERE ) result ( Node_Kind )
  ! Return the kind of the tree node at WHERE
    integer, intent(in) :: WHERE
    node_kind = the_tree(where) % kind
  end function NODE_KIND_I

  pure integer function NODE_KIND_TX ( WHERE ) result ( Node_Kind )
  ! Return the kind ofthe tree node at WHERE
    type(tx), intent(in) :: WHERE
    node_kind = the_tree(where%i) % kind
  end function NODE_KIND_TX

  pure integer function NSONS_I ( WHERE ) result ( Nsons )
  ! Return the node id of the tree node at WHERE
    integer, intent(in) :: WHERE
    nsons = the_tree(where) % nsons
  end function NSONS_I

  pure integer function NSONS_TX ( WHERE ) result ( Nsons )
  ! Return the node id of the tree node at WHERE
    type(tx), intent(in) :: WHERE
    nsons = the_tree(where%i) % nsons
  end function NSONS_TX

  subroutine POP ( N )
  ! Delete the top N nodes from the orchard stack
    integer, intent(in) :: N
    tree_sp = tree_sp + n
    n_tree_stack = n_tree_stack - n
  end subroutine POP

  recursive subroutine PRINT_SUBTREE_I ( SUBROOT, DEPTH, DUMP_DECOR, TYPE_NAME )
  ! Print the subtree rooted at SUBROOT, starting with DEPTH leading
  ! dots.  Display the decoration of each tree node if DUMP_DECOR is
  ! present and .true.
    integer, intent(in) :: SUBROOT
    integer, intent(in) :: DEPTH
    logical, intent(in), optional :: DUMP_DECOR
    optional :: Type_Name
    interface
      integer function Type_Name ( Decor )
      ! Return the string index to print for the decoration
        integer, intent(in) :: Decor ! Tree(where)%Decor
      end function Type_Name
    end interface
    integer :: I, MyRoot
    myRoot = subroot
    if ( myRoot < 0 ) myRoot = tree_sp + 1
    call output ( myRoot, 5 )
    call output ( ':' )
    call dump_tree_node ( myRoot, depth, type_name=type_name )
    if ( present(dump_decor) ) then
      if ( dump_decor ) then
      ! In Fortran 2000:
      ! write ( prunit, '(" (", i0, ") ', advance="no") decoration(myRoot)
        call output ( ' (' )
        call output ( decoration(myRoot) )
        call output ( ') ')
      end if
    end if
    call newLine
    if ( the_tree(myRoot)%kind == internal ) then
      do i = 1, nsons(myRoot)
        call print_subtree ( subtree(i,myRoot), depth+1, dump_decor, &
          & type_name=type_name )
      end do
    end if
  end subroutine PRINT_SUBTREE_I

  recursive subroutine PRINT_SUBTREE_TX ( SUBROOT, DEPTH, DUMP_DECOR, TYPE_NAME )
  ! Print the subtree rooted at SUBROOT, starting with DEPTH leading
  ! dots.  Display the decoration of each tree node if DUMP_DECOR is
  ! present and .true.
    type(tx), intent(in) :: SUBROOT
    integer, intent(in) :: DEPTH
    logical, intent(in), optional :: DUMP_DECOR
    interface
      integer function Type_Name ( Decor )
      ! Return the string index to print for the decoration
        integer, intent(in) :: Decor ! Tree(where)%Decor
      end function Type_Name
    end interface
    call print_subtree ( subroot%i, depth, dump_decor, type_name )
  end subroutine PRINT_SUBTREE_TX

  subroutine PUSH_PSEUDO_TERMINAL_INTEGER ( SUB_ROSA, SOURCE, DECOR, FILE, &
    & CLASS )
  ! Push the pseudo-terminal with string index SUB_ROSA, source SOURCE
  ! and decoration DECOR onto the tree stack.  If DECOR is absent, use
  ! null_tree.
    integer, intent(in) :: SUB_ROSA
    integer, intent(in) :: SOURCE
    integer, intent(in), optional :: DECOR
    integer, intent(in), optional :: FILE
    integer, intent(in), optional :: CLASS ! of string
    integer :: MY_CLASS, MY_DECOR, MY_FILE
    my_decor = null_tree
    if ( present(decor) ) my_decor = decor
    my_file = 0
    if ( present(file) ) my_file = file
    my_class = symbol(sub_rosa)
    if ( present(class) ) my_class = class

    if ( tree_sp <= tree_point ) then
      call tree_error ( no_tree_space, null_tree )
!     call double_tree
    end if
                                 ! node                   decor
    the_tree(tree_sp) = tree_node( tree_map(my_class), my_decor, &
                             ! file nsons  source  kind    sub_rosa
                               my_file, 0, source, pseudo, sub_rosa )
    tree_sp = tree_sp - 1     ! push tree stack
    n_tree_stack = n_tree_stack + 1
  end subroutine PUSH_PSEUDO_TERMINAL_INTEGER

  subroutine PUSH_PSEUDO_TERMINAL_WHERE (SUB_ROSA, WHERE, DECOR, &
    & CLASS )
  ! Push the pseudo-terminal with string index SUB_ROSA, source SOURCE
  ! and decoration DECOR onto the tree stack.  If DECOR is absent, use
  ! null_tree.
    use Lexer_Core, only: Where_t
    integer, intent(in) :: SUB_ROSA
    type(where_t), intent(in) :: WHERE
    integer, intent(in), optional :: DECOR
    integer, intent(in), optional :: CLASS ! of string
    call push_pseudo_terminal ( sub_rosa, where%source, decor, where%file, class )
  end subroutine PUSH_PSEUDO_TERMINAL_WHERE

  subroutine REPLACE_SONS_I ( T, K, M, U )
  ! Replace sons k .. m of t with the tree node at u.  This will leave an
  ! empty space if k < m.
    integer, intent(in) :: T            ! The tree node whose sons are to
                                        ! be replaced
    integer, intent(in) :: K            ! The first son of T to be replaced
    integer, intent(in) :: M            ! The last son of T to be replaced
    integer, intent(in) :: U            ! The new tree node
    integer :: I, J                     ! Tree nodes to copy to and from
    integer :: MK, MM                   ! my K, my M
    integer :: N, NN                    ! NSONS(T), new value for NSONS(T)
    i = subtree(k,t)
    the_tree(i) = the_tree(u)
    mk = k
    mm = m
    if ( mm < mk ) then
      mm = mk
    else if ( mk < m ) then
      n = nsons(t)
      nn = n - ( mm-mk )
      do while ( mm < n )
        mk = mk + 1
        mm = mm + 1
        i = subtree( mk, t )
        j = subtree( mm, t )
        the_tree(i) = the_tree(j)
      end do
      mk = mk + 1
      i = subtree( mk, t )
      do while ( mk <= n )
      ! we can't just use "i = subtree(k,t)" because that
      !wouldn't convert the "more" nodes to "empty" nodes. *)
        if ( the_tree(i)%kind == more ) then
          the_tree(i)%kind = empty
          i = the_tree(i)%link
        end if
        the_tree(i)%kind = empty
        i = i + 1
        mk = mk + 1
      end do
      call set_nsons( t, nn )
      ! decrease treept if we make empties at the end of the tree
      do while ( the_tree(tree_sp)%kind == empty )
        tree_sp = tree_sp - 1
      end do
    end if
  end subroutine REPLACE_SONS_I

  subroutine REPLACE_SONS_TX ( T, K, M, U )
  ! Replace sons k .. m of t with the tree node at u.  This will leave an
  ! empty space if k < m.
    type(tx), intent(in) :: T           ! The tree node whose sons are to
                                        ! be replaced
    integer, intent(in) :: K            ! The first son of T to be replaced
    integer, intent(in) :: M            ! The last son of T to be replaced
    type(tx), intent(in) :: U           ! The new tree node
    call replace_sons ( t%i, k, m, u%i )
  end subroutine REPLACE_SONS_TX

  pure integer function SOURCE_REF_I ( WHERE ) result ( Source_Ref )
  ! Return the SOURCE field of the tree node at WHERE
    integer, intent(in) :: WHERE
    source_ref = the_tree(where) % source
  end function SOURCE_REF_I

  pure integer function SOURCE_REF_TX ( WHERE ) result ( Source_Ref )
  ! Return the SOURCE field of the tree node at WHERE
    type(tx), intent(in) :: WHERE
    source_ref = the_tree(where%i) % source
  end function SOURCE_REF_TX

  pure integer function STACK_FILE ( ) result ( The_File )
  ! Return the SOURCE field of the top stack frame
    the_file = the_tree(tree_sp+1) % file
  end function STACK_FILE

  pure integer function STACK_SOURCE_REF ( ) result ( Source_Ref )
  ! Return the SOURCE field of the top stack frame
    source_ref = the_tree(tree_sp+1) % source
  end function STACK_SOURCE_REF

  integer function STACK_SUBTREE ( WHICH )
  ! Return the root of the WHICH'th subtree of the node atop the stack
    integer, intent(in) :: WHICH
    stack_subtree = subtree(which, tree_sp+1)
  end function STACK_SUBTREE

  type(tx) function STACK_SUBTREE_TX ( WHICH )
  ! Return the root of the WHICH'th subtree of the node atop the stack
    integer, intent(in) :: WHICH
    stack_subtree_tx%i = subtree(which, tree_sp+1)
  end function STACK_SUBTREE_TX

  integer function STACK_SUB_ROSA ( )
  ! Return the sub_rosa index of the top stack frame
    stack_sub_rosa = the_tree(tree_sp+1) % link
  end function STACK_SUB_ROSA

  integer function SUB_ROSA_I ( WHERE ) result ( Sub_Rosa )
  ! Return the sub_rosa string index from the tree node at WHERE
    integer, intent(in) :: WHERE
    if ( the_tree(where) % kind /= pseudo ) then
      call tree_error ( not_pseudo, where )
    end if
    sub_rosa = the_tree(where) % link
  end function SUB_ROSA_I

  integer function SUB_ROSA_TX ( WHERE )
  ! Return the sub_rosa string index from the tree node at WHERE
    type(tx), intent(in) :: WHERE
    sub_rosa_tx = sub_rosa ( where % i )  
  end function SUB_ROSA_TX

  integer function SUBTREE_I ( WHICH, WHERE ) result ( Subtree )
  ! Return the root of the WHICH'th subtree of the tree node at WHERE
  ! The zero'th subtree of WHERE is WHERE
    integer, intent(in) :: WHICH
    integer, intent(in) :: WHERE
    if ( the_tree(where) % kind == pseudo ) then
      call tree_error ( is_pseudo, where )
    end if
    if ( which == 0 ) then
      subtree = where
      return
    end if
    if ( which < 0 .or. which > the_tree(where) % nsons ) then
      call tree_error ( no_more_sons, where, which, the_tree(where) % nsons )
    end if
    subtree = the_tree(where) % link + which - 1
  end function SUBTREE_I

  type(tx) function SUBTREE_TX ( WHICH, WHERE )
  ! Return the root of the WHICH'th subtree of the tree node at WHERE
  ! The zero'th subtree of WHERE is WHERE
    integer, intent(in) :: WHICH
    type(tx), intent(in) :: WHERE
    subtree_tx%i = Subtree ( which, where%i )
  end function SUBTREE_TX

  pure integer function THE_FILE_I ( WHERE ) result ( The_File )
  ! Return the FILE field of the tree node at WHERE
    integer, intent(in) :: WHERE
    the_file = the_tree(where) % file
  end function THE_FILE_I

  pure integer function THE_FILE_TX ( WHERE ) result ( The_File )
  ! Return the FILE field of the tree node at WHERE
    type(tx), intent(in) :: WHERE
    the_file = the_tree(where%i) % file
  end function THE_FILE_TX

  pure integer function TREE_NODE_NAME_I ( WHERE )
  ! Get the string index of the name of the tree node at WHERE
    integer, intent(in) :: WHERE
    tree_node_name_i = tree_texts(the_tree(where) % node)
  end function TREE_NODE_NAME_I

  pure integer function TREE_NODE_NAME_TX ( WHERE )
  ! Get the string index of the name of the tree node at WHERE
    type(tx), intent(in) :: WHERE
    tree_node_name_tx = tree_texts(the_tree(where%i) % node)
  end function TREE_NODE_NAME_TX

  pure integer function TREE_TEXT ( TREE_NODE )
  ! Return the string index of the text of the tree node id TREE_NODE
    integer, intent(in) :: TREE_NODE
    tree_text = tree_texts ( tree_node )
  end function TREE_TEXT

  function Where_I ( Tree ) result ( Where )
  ! Return the Where_T structure of the text of the tree node at WHERE
    use Lexer_Core, only: Where_T
    integer, intent(in) :: Tree ! Tree node index
    type(where_t) :: Where
    where%source = the_tree(tree)%source
    where%file = the_tree(tree)%file
  end function Where_I

  function Where_TX ( Tree ) result ( Where )
  ! Return the Where_T structure of the text of the tree node at WHERE
    use Lexer_Core, only: Where_T
    type(tx), intent(in) :: Tree ! Tree node index
    type(where_t) :: Where
    where%source = the_tree(tree%i)%source
    where%file = the_tree(tree%i)%file
  end function Where_TX

! =====     Private procedures     =======================================

  subroutine DOUBLE_TREE ( STATUS )
  ! Double the tree space
    integer, intent(out), optional :: STATUS
    type(tree_node), allocatable :: Old_Tree(:)
    integer :: New_SP         ! New tree stack pointer
    integer :: STAT           ! From ALLOCATE

    allocate ( old_tree(size(the_tree)), stat=stat )
    if ( stat /= 0 ) then
      if ( present(status) ) then
        status = stat
        return
      end if
      call io_error ( 'TREE%DOUBLE_TREE-E- Unable to allocate tree', &
        stat, ' ' )
      call StopWithErrorMsg ( 'TREE%DOUBLE_TREE-E- Unable to allocate old tree' )
    end if
    old_tree = the_tree
    deallocate ( the_tree )
    allocate ( the_tree(2*size(old_tree)), stat=stat )
    if ( stat /= 0 ) then
      if ( present(status) ) then
        status = stat
        return
      end if
      call io_error ( 'TREE%DOUBLE_TREE-E- Unable to allocate tree', &
        stat, ' ' )
      call StopWithErrorMsg ( 'TREE%DOUBLE_TREE-E- Unable to allocate new tree' )
    end if
    new_sp = size(the_tree) - ( size(old_tree) - tree_sp )
    the_tree(:tree_point) = old_tree(:tree_point)
    the_tree(new_sp:) = old_tree(tree_sp:)
    deallocate( old_tree )
    tree_sp = new_sp
  end subroutine DOUBLE_TREE

  subroutine POP_TO_TREE ( NSONS )
  ! Pop the top nSons nodes from the tree stack into the end of the tree
  ! space.
    integer, intent(in) :: NSONS        ! How many sons to copy
    integer :: I, N_COPY, N_EXCH
    type(tree_node) :: TEMP_NODE
    ! Copy tree stack to tree.  Reverse the order.  Exchange the part that
    ! overlaps, if any.
    if ( tree_point + nsons <= tree_sp ) then; n_copy = nsons
    else; n_copy = tree_sp - tree_point; end if
    n_exch = n_copy + ( nsons - n_copy ) / 2
    do i = 1, n_copy ! We don't use an array section because it would
                     ! create an array temp
      the_tree(tree_point+i) = the_tree(tree_sp + nsons - i + 1)
    end do
    do i = n_copy+1, n_exch
      temp_node = the_tree(tree_point + i)
      the_tree(tree_point+i) = the_tree(tree_sp + nsons - i + 1)
      the_tree(tree_sp + nsons - i + 1) = temp_node
    end do
    tree_sp = tree_sp + nsons ! Stack pop
    tree_point = tree_point + nsons
  end subroutine POP_TO_TREE

  subroutine SET_NSONS ( T, NSONS )
  ! Set the number of sons.  This isn't straight-forward in the case that
  ! the last son is a "more" node.
    integer, intent(in) :: T            ! The tree node to change
    integer, intent(in) :: NSONS        ! The new number of sons
    integer :: N, NS, S, TT
    ns = nsons
    tt = t
    do
      n = the_tree(tt)%nsons
      if ( ns < n ) then
        the_tree(tt)%nsons = ns
    exit
      end if
      s = the_tree(tt)%link + n - 1
      if ( the_tree(s)%kind /= more ) then
        ! we better have increased nsons by storing adjacent to old ones!
        the_tree(t)%nsons = n
    exit
      end if
      tt = s
      ns = ns - ( n-1 )
    end do
  end subroutine SET_NSONS

  subroutine TREE_ERROR ( WHY, WHERE, Which, How_Many )
  ! Print a message WHY there's an error at WHERE in the tree.
    integer, intent(in) :: WHY
    integer, intent(in) :: WHERE
    integer, intent(in), optional :: Which
    integer, intent(in), optional :: How_Many
    if ( where == null_tree ) then
      call error_intro ( compiler, 0 )
    else
      call error_intro ( compiler, source_ref(where) )
    end if
    if ( why == is_pseudo .or. why == no_more_sons .or. why == not_pseudo ) &
      & call output ( where, before="Tree node at ")
    select case ( why )
    case ( is_pseudo );     call output ( " is pseudo-terminal.", &
                                          advance="yes" )
    case ( no_more_sons);
      call output ( " has too few sons" )
      if ( present(which) ) call output ( which, before=".  Asking for subtree " )
      if ( present(how_Many) ) call output ( how_Many, before=" out of " )
      call output ( '.', advance="yes" )
    case ( not_pseudo );    call output ( " is not pseudo-terminal.", &
                                          advance="yes" )
    case ( no_tree_space )
      call output ( "No Tree space remains.", advance="yes" )
      call output ( "Unfortunately, 'tree' doesn't know how to increase its space", &
        & advance="yes" )
    case ( underflow );     call output ( "Tree stack underflow.", &
                                               advance="yes" )
    end select
    call StopWithErrorMsg ( 'Possibly an error in the source code' )
  end subroutine TREE_ERROR

  ! ------------------ Private ---------------------
  ! ------------ StopWithErrorMsg ------------
  subroutine StopWithErrorMsg ( Message )
    ! Print Message, dump calling stack (if any) and stop
    character (len=*), intent(in) :: Message ! Line of text
    ! Executable
    call PrintItOut( message, MLSMSG_ERROR, exitStatus = 1  )
  end subroutine StopWithErrorMsg

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module TREE

! $Log$
! Revision 2.33  2018/08/03 20:25:29  vsnyder
! Repair two comments
!
! Revision 2.32  2015/09/17 22:47:29  pwagner
! Added Node_in_tree
!
! Revision 2.31  2014/05/30 02:43:38  vsnyder
! Improve no_more_sons error message
!
! Revision 2.30  2014/05/20 22:16:57  vsnyder
! More functions to access the stack
!
! Revision 2.29  2014/02/21 19:21:31  vsnyder
! Use First_Tree_Node instead of N_Eof as the first tree node
!
! Revision 2.28  2014/01/11 01:41:02  vsnyder
! Decruftification
!
! Revision 2.27  2014/01/08 21:07:58  vsnyder
! Get_Tree_Node_Name_TX ought to get the name, not display it
!
! Revision 2.26  2013/12/12 02:01:17  vsnyder
! Add 'type_name' dummy procedure to Dump_Tree_Node
!
! Revision 2.25  2013/11/26 22:46:52  vsnyder
! Add Dump_Stack, class of string in Push_Pseudo_Terminal, Stack_Sub_Rosa
!
! Revision 2.24  2013/10/09 01:09:42  vsnyder
! Add some routines, spiff up dumps
!
! Revision 2.23  2013/10/02 01:31:06  vsnyder
! Add Dump_Top_Stack, Dump_Top_Stack_Name
!
! Revision 2.22  2013/09/30 23:59:06  vsnyder
! Default initializer for tx%i, add decoration_tx_tx
!
! Revision 2.21  2013/09/30 23:03:04  vsnyder
! Add TX type for tree index and generics to use it
!
! Revision 2.20  2013/09/24 23:09:15  vsnyder
! Replace Source with Where_t, add Where function
!
! Revision 2.19  2013/09/19 23:25:57  vsnyder
! Add some tracing
!
! Revision 2.18  2013/09/12 03:12:46  vsnyder
! Add Advance and Before to Dump_Tree_Node_Name
!
! Revision 2.17  2013/08/28 00:36:21  pwagner
! Moved more stuff from MLSMessage down to PrintIt module
!
! Revision 2.16  2012/03/15 22:46:37  vsnyder
! Make nsons pure, some cannonball polishing
!
! Revision 2.15  2011/04/18 19:26:11  vsnyder
! Add Dump_Tree_Node_Name
!
! Revision 2.14  2010/04/28 00:13:47  pwagner
! Replaced bare stop with StopWithErrorMsg
!
! Revision 2.13  2009/09/29 23:22:44  vsnyder
! Arguments reversed in call to error_intro
!
! Revision 2.12  2009/06/23 18:25:44  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.11  2008/09/04 00:45:51  vsnyder
! Add PRINT to not_used_here
!
! Revision 2.10  2006/07/28 02:00:17  vsnyder
! Pure cannonball polishing
!
! Revision 2.9  2006/02/23 00:56:43  vsnyder
! Add source line and column to tree node dump
!
! Revision 2.8  2005/06/22 17:25:51  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.7  2004/05/28 23:14:35  vsnyder
! Allow print_subtree to start at the top of the stack with a negative 'root'
!
! Revision 2.6  2003/05/12 20:54:16  vsnyder
! Correct a subtle bug in SUBTREE that probably doesn't affect us
!
! Revision 2.5  2002/10/08 00:09:15  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2001/04/05 00:55:27  vsnyder
! Try to write 'increase tree size automatically' code -- failed -- try again later
!
! Revision 2.3  2001/02/23 00:33:11  vsnyder
! Move source_ref to the beginning of the subtree (from the end)
!
! Revision 2.2  2001/02/07 18:40:22  vsnyder
! Add a "decor" argument to "build_tree".
!
! Revision 2.1  2000/10/11 18:33:25  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:51  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
