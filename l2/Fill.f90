! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Fill                     ! Create vectors and fill them.
!=============================================================================

  use GriddedData, only: GriddedData_T
  use INIT_TABLES_MODULE, only: S_VECTOR
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: L1BInfo_T
  use OUTPUT_M, only: OUTPUT
  use QuantityTemplates, only: QuantityTemplate_T
  use TOGGLES, only: GEN, LEVELS, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, DUMP_TREE_NODE, NODE_ID, NSONS, &
    & SOURCE_REF, SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_NAMED
  use VectorsModule, only: AddVectorToDatabase, CreateVector, Dump, Vector_T, &
    & VectorTemplate_T

  implicit none
  private
  public :: MLSL2Fill
  
  ! -----     Private declarations     ---------------------------------

  integer, private :: ERROR

  ! Error codes for "announce_error"
  integer, parameter :: WRONG_NUMBER = 1     ! of fields of a VECTOR command

  !---------------------------- RCS Ident Info -------------------------------
  character (len=256) :: Id = &
       "$id: fill.f90,v 1.1 2000/01/21 21:04:06 livesey Exp $"
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module performs the Fill operation in the Level 2 software.  
  ! This takes a vector template, and creates and fills an appropriate vector

contains ! =====     Public Procedures     =============================

  !---------------------------------------------------  MLSL2Fill  -----

  subroutine MLSL2Fill ( root, l1bInfo, aprioriData, vectorTemplates, vectors, &
    & qtyTemplates )

  ! This is the main routine for the module.  It parses the relevant lines
  ! of the l2cf and works out what to do.

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the FILL section in the AST
    type (L1BInfo_T), intent(in) :: l1bInfo
    type (GriddedData_T), dimension(:), pointer :: aprioriData
    type (VectorTemplate_T), dimension(:), pointer :: vectorTemplates
    type (Vector_T), dimension(:), pointer :: vectors
    type (QuantityTemplate_T), dimension(:), pointer :: qtyTemplates

    ! Local variables
    integer :: I, J                ! Loop indices for section, spec
    integer :: KEY                 ! Definitely n_named
    type (Vector_T) :: newVector
    integer :: SON                 ! Of root, an n_spec_args or a n_named
    integer :: templateIndex       ! In the template database
    integer :: vectorIndex         ! In the vector database
    integer :: vectorName          ! Sub-rosa index

    ! Executable code

    if ( toggle(gen) ) call trace_begin ( "MLSL2Fill", root )

    error = 0
    templateIndex = -1
    vectorIndex = -1

    ! Loop over the lines in the configuration file

    do i = 2, nsons(root)-1 ! Skip the section name at begin and end
      son = subtree(i,root)
      if ( node_id(son) == n_named ) then ! Is spec labeled?
        key = subtree(2,son)
        vectorName = sub_rosa(subtree(1,son))
      else
        key = son
        vectorName = 0
      end if

      ! Node_id(key) is now n_spec_args.

      select case( decoration(subtree(1,decoration(subtree(1,key)))) )
      case ( s_vector )
        if ( nsons(key) /= 2 ) call announce_error ( son, wrong_number )
        templateIndex = decoration(decoration(subtree(2,subtree(2,key))))

        ! Create the vector, and add it to the database.

        call decorate ( key, AddVectorToDatabase ( vectors, &
          & CreateVector ( vectorName, vectorTemplates(templateIndex), &
            & qtyTemplates ) ) )

        ! That's the end of the create operation

!     case ( s_fill )

      case default ! Can't get here if tree_checker worked correctly
      end select
    end do

    if ( toggle(gen) ) then
      if ( levels(gen) > 0 ) then
        call dump ( vectors )
      end if
      call trace_end ( "MLSL2Fill" )
    end if
  end subroutine MLSL2Fill

! =====     Private Procedures     =====================================

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref(where) )
    call output ( ': ' )
    select case ( code )
    case ( wrong_number )
      call output ( "The " );
      call dump_tree_node ( where, 0 )
      call output ( " command does not have exactly one field.", advance='yes' )
    end select
  end subroutine ANNOUNCE_ERROR

!=============================================================================
end module Fill
!=============================================================================

!
! $Log$
! Revision 2.2  2000/09/11 19:52:51  ahanzel
! Removed old log entries in file.
!
! Revision 2.1  2000/09/08 22:55:56  vsnyder
! Revised to use the tree output by the parser
!
!
