! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module TREE

  use ERROR_HANDLER, only: COMPILER, ERROR_INTRO
  use MACHINE, only: IO_ERROR
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: DISPLAY_STRING, LOOKUP_AND_INSERT
  use SYMBOL_TABLE, only: SET_SYMBOL, SYMBOL
  use SYMBOL_TYPES, only: T_NULL
  use TOGGLES, only: CON, TOGGLE
  use TREE_TYPES, only: LAST_TREE_NODE, N_EOF, N_NULL, TREE_INIT, TREE_MAP
  implicit NONE
  private

  public :: ADD_SONS_FROM_STACK, ALLOCATE_TREE, BUILD_TREE, COPY_TO_STACK
  public :: DEALLOCATE_TREE, DECORATE, DECORATION, DELETE_TREE
  public :: DELETE_TREE_STACK, DUMP_TREE_NODE, INIT_TREE, INSERT_NODE
  public :: NODE_ID, NODE_KIND, NSONS, POP, PRINT_SUBTREE
  public :: PUSH_PSEUDO_TERMINAL, REPLACE_SONS, SOURCE_REF, STACK_SUBTREE
  public :: SUB_ROSA, SUBTREE, TREE_TEXT

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

  type :: TREE_NODE
    integer :: NODE      ! What kind of tree node
    integer :: DECOR     ! Decoration, an integer or tree node index
    integer :: NSONS     ! How many sons
    integer :: SOURCE    ! 256*line + column
    integer :: KIND      ! PSEUDO, INTERNAL, MORE or OTHER
    integer :: LINK      ! Sub_rosa if PSEUDO, Left_son if INTERNAL, next
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

  integer, save, private :: TREE_TEXTS(N_EOF:LAST_TREE_NODE)

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine ADD_SONS_FROM_STACK ( T, N )
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
  end subroutine ADD_SONS_FROM_STACK

  subroutine ALLOCATE_TREE ( N_TREE, STATUS )
  ! Allocate a TREE array with N_TREE elements
    integer, intent(in) :: N_TREE
    integer, intent(out), optional :: STATUS ! from ALLOCATE
    integer :: STAT
    if ( allocated(the_tree) ) then; deallocate ( the_tree ); end if
    allocate ( the_tree(0:n_tree), stat=stat )
    if ( stat /= 0 ) then
      if ( present(status) ) then
        status = stat
        return
      end if
      call io_error ( 'TREE%ALLOCATE_TREE-E- Unable to allocate tree', &
        stat, ' ' )
      stop
    end if
    call init_tree
  end subroutine ALLOCATE_TREE

  subroutine BUILD_TREE ( NEW_NODE, NSONS, DECORATION )
  ! Pop NSONS nodes from the orchard stack to the tree.  Build a new tree
  ! node of type NEW_NODE over them.  Leave the new node on the stack at
  ! TREE_SP+1.
    integer, intent(in) :: NEW_NODE
    integer, intent(in) :: NSONS
    integer, intent(in), optional :: DECORATION
    integer :: LEFT_SON, MY_DECOR, SOURCE_REF
    my_decor = null_tree
    if ( present(decoration) ) my_decor = decoration
    if ( tree_sp + nsons > ubound(the_tree,1) ) then
      call tree_error ( underflow, null_tree )
    end if
    if ( nsons == 0 ) then; source_ref = 0
    else; source_ref = the_tree(tree_sp+nsons)%source; end if
    left_son = tree_point + 1
    call pop_to_tree ( nsons )
    ! Push new node onto the stack
    if ( tree_sp <= tree_point ) then
      call tree_error ( no_tree_space, null_tree )
!     call double_tree
    end if
                                   ! node      decor   nsons  source
    the_tree(tree_sp) = tree_node( new_node, my_decor, nsons, source_ref, &
                                 !    kind   left_son
                                   internal, left_son )
    tree_sp = tree_sp - 1     ! Stack push
    n_tree_stack = n_tree_stack + 1 - nsons
    return
  end subroutine BUILD_TREE

  subroutine COPY_TO_STACK ( T, K, M )
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
    return
  end  subroutine COPY_TO_STACK

  subroutine DEALLOCATE_TREE
    deallocate( the_tree )
    return
  end subroutine DEALLOCATE_TREE

  subroutine DECORATE ( WHERE, DECORATION )
  ! Decorate tree WHERE with DECORATION
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
  end subroutine DECORATE

  integer function DECORATION ( WHERE )
  ! Return the decoration of the tree WHERE
    integer, intent(in) :: WHERE
    decoration = the_tree(where) % decor
  end function DECORATION

  subroutine DELETE_TREE ( WHERE )
  ! Discard everything at WHERE and above in the tree
    integer, intent(in) :: WHERE
    tree_point = where - 1
  end subroutine DELETE_TREE

  subroutine DELETE_TREE_STACK
  ! Discard everything in the tree stack
    tree_sp = ubound(the_tree,1)
    n_tree_stack = 0
  end subroutine DELETE_TREE_STACK

  subroutine DUMP_TREE_NODE ( WHERE, INDENT, ADVANCE )
  ! Indent INDENT spaces, then dump the tree WHERE
    integer, intent(in) :: WHERE
    integer, intent(in) :: INDENT
    character(len=*), intent(in), optional :: ADVANCE
    integer :: I
    do i = 1, indent; call output ( '.' ); end do
    select case ( the_tree(where) % kind )
    case ( empty )
      call output ( 'empty', advance )
    case ( more )
      call output ( 'more, nsons = ' )
      call output ( the_tree(where) % nsons )
      call output ( ', next = ' )
      call output ( the_tree(where) % link, advance=advance )
    case ( pseudo, internal )
      if ( the_tree(where) % kind == pseudo ) then
        call display_string ( tree_texts(the_tree(where) % node) )
        call output ( ' ' )
        call display_string ( sub_rosa(where) )
      else
        call display_string ( tree_texts(the_tree(where) % node) )
      end if
      if ( the_tree(where)%decor /= null_tree ) then
        call output ( ' decor=' )
        call output ( the_tree(where) % decor )
      end if
      call output ( '', advance )
    end select
    return
  end subroutine DUMP_TREE_NODE

  subroutine INIT_TREE
    logical :: FOUND     ! Did lookup_and_insert find it?
    integer :: I         ! Loop inductor
    integer :: WHERE     ! Where did lookup_and_insert find it?
    do i = n_eof, last_tree_node
      call tree_init (i)
      call lookup_and_insert ( where, found, .false. )
      ! It's OK if it found one -- maybe it's a terminal text, too.
      if ( .not. found ) then; call set_symbol(where, t_null); end if
      tree_texts(i) = where
    end do
    tree_point = null_tree
    call delete_tree_stack
                                      ! node  decor  nsons  source kind
    the_tree(tree_point) = tree_node( n_null, null_tree, 0, 0, internal, &
                                      ! left_son
                                      null_tree )
    return
  end subroutine INIT_TREE

  subroutine INSERT_NODE ( NEW_NODE, T, K, M )
  ! Insert a node having ID = newNode as the k'th son of t.  Make the sons
  ! k..m of t sons of newNode.  Reduce the number of sons of t to
  ! n-k+m-1.  Don't use this procedure during parsing:  It checks nTree,
  ! not treesp!
    integer, intent(in) :: NEW_NODE     ! ID of the node to be inserted
    integer, intent(in) :: T            ! The parent of the inserted node
    integer, intent(in) :: K            ! Which sone of T NEW_NODE is to be
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
    return
  end subroutine INSERT_NODE

  integer function NODE_ID ( WHERE )
  ! Return the node id of the tree node at WHERE
    integer, intent(in) :: WHERE
    node_id = the_tree(where) % node
    return
  end function NODE_ID

  integer function NODE_KIND ( WHERE )
  ! Return the kind of tree(where).
    integer, intent(in) :: WHERE
    node_kind = the_tree(where) % kind
    return
  end function NODE_KIND

  integer function NSONS ( WHERE )
  ! Return the node id of the tree node at WHERE
    integer, intent(in) :: WHERE
    nsons = the_tree(where) % nsons
    return
  end function NSONS

  subroutine POP ( N )
  ! Delete the top N nodes from the orchard stack
    integer, intent(in) :: N
    tree_sp = tree_sp + n
    n_tree_stack = n_tree_stack - n
    return
  end subroutine POP

  recursive subroutine PRINT_SUBTREE ( SUBROOT, DEPTH, DUMP_DECOR )
  ! Print the subtree rooted at SUBROOT, starting with DEPTH leading
  ! dots.  Display the decoration of each tree node if DUMP_DECOR is
  ! present and .true.
    integer, intent(in) :: SUBROOT
    integer, intent(in) :: DEPTH
    logical, intent(in), optional :: DUMP_DECOR
    integer :: I
    call output ( subroot, 5 )
    call output ( ':' )
    call dump_tree_node ( subroot, depth )
    if ( present(dump_decor) ) then
      if ( dump_decor ) then
      ! In Fortran 2000:
      ! write ( prunit, '(" (", i0, ") ', advance="no") decoration(subroot)
        call output ( ' (' )
        call output ( decoration(subroot) )
        call output ( ') ')
      end if
    end if
    call output ( '', advance='yes' )
    if ( the_tree(subroot)%kind == internal ) then
      do i = 1, nsons(subroot)
        call print_subtree ( subtree(i,subroot), depth+1, dump_decor )
      end do
    end if
    return
  end subroutine PRINT_SUBTREE

  subroutine PUSH_PSEUDO_TERMINAL ( SUB_ROSA, SOURCE, DECOR )
  ! Push the pseudo-terminal with string index SUB_ROSA, source SOURCE
  ! and decoration DECOR onto the tree stack.  If DECOR is absent, use
  ! null_tree.
    integer, intent(in) :: SUB_ROSA
    integer, intent(in) :: SOURCE
    integer, intent(in), optional :: DECOR
    integer :: MY_DECOR
    my_decor = null_tree
    if ( present(decor) ) my_decor = decor
    if ( tree_sp <= tree_point ) then
      call tree_error ( no_tree_space, null_tree )
!     call double_tree
    end if
                                   ! node                       decor
    the_tree(tree_sp) = tree_node( tree_map(symbol(sub_rosa)), my_decor, &
                             ! nsons  source  kind    sub_rosa
                                   0, source, pseudo, sub_rosa )
    tree_sp = tree_sp - 1     ! push tree stack
    n_tree_stack = n_tree_stack + 1
    return
  end subroutine PUSH_PSEUDO_TERMINAL

  subroutine REPLACE_SONS ( T, K, M, U )
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
    return
  end subroutine REPLACE_SONS

  integer function SOURCE_REF ( WHERE )
  ! Return the SOURCE field of the tree node at WHERE
    integer, intent(in) :: WHERE
    source_ref = the_tree(where) % source
    return
  end function SOURCE_REF

  integer function STACK_SUBTREE ( WHICH )
  ! Return the root of the WHICH'th subtree of the node atop the stack
    integer, intent(in) :: WHICH
    stack_subtree = subtree(which, tree_sp+1)
    return
  end function STACK_SUBTREE

  integer function SUB_ROSA ( WHERE )
  ! Return the sub_rosa string pointer from the tree node at WHERE
    integer, intent(in) :: WHERE
    if ( the_tree(where) % kind /= pseudo ) then
      call tree_error ( not_pseudo, where )
    end if
    sub_rosa = the_tree(where) % link
    return
  end function SUB_ROSA

  integer function SUBTREE ( WHICH, WHERE )
  ! Return the root of the WHICH'th subtree of the tree node at WHERE
  ! The zero'th subtree of WHERE is WHERE
    integer, intent(in) :: WHICH
    integer, intent(in) :: WHERE
    if ( the_tree(where) % kind == pseudo ) then
      call tree_error ( is_pseudo, where )
    end if
    if ( which == 0 ) then
      subtree = which
      return
    end if
    if ( which < 0 .or. which > the_tree(where) % nsons ) then
      call tree_error ( no_more_sons, where )
    end if
    subtree = the_tree(where) % link + which - 1
    return
  end function SUBTREE

  integer function TREE_TEXT ( TREE_NODE )
  ! Return the string index of the text of the tree node
    integer, intent(in) :: TREE_NODE
    tree_text = tree_texts(tree_node)
  end function TREE_TEXT
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
      stop
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
      stop
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
    return
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
    return
  end subroutine SET_NSONS

  subroutine TREE_ERROR ( WHY, WHERE )
  ! Print a message WHY there's an error at WHERE in the tree.
    integer, intent(in) :: WHY
    integer, intent(in) :: WHERE
    if ( where == null_tree ) then
      call error_intro ( 0, compiler )
    else
      call error_intro ( source_ref(where), compiler )
    end if
    if ( why == is_pseudo .or. why == no_more_sons .or. why == not_pseudo ) &
    then
      ! write ( prunit, '(a,i0)', advance="no" ) where ! in Fortran 2000
      call output ( "Tree node at ")
      call output ( where )
    end if
    select case ( why )
    case ( is_pseudo );     call output ( " is pseudo-terminal.", &
                                          advance="yes" )
    case ( no_more_sons);   call output ( " has too few sons.", &
                                          advance="yes" )
    case ( not_pseudo );    call output ( " is not pseudo-terminal.", &
                                          advance="yes" )
    case ( no_tree_space )
      call output ( "No Tree space remains.", advance="yes" )
      call output ( "Unfortunately, 'tree' doesn't know how to increase its space", &
        & advance="yes" )
    case ( underflow );     call output ( "Tree stack underflow.", &
                                               advance="yes" )
    end select
    stop
  end subroutine TREE_ERROR

end module TREE

! $Log$
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
