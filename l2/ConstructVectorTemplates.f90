! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ConstructVectorTemplates ! Construct a template for a vector
!=============================================================================

  use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
  use INIT_TABLES_MODULE, only: F_QUANTITIES, F_SIGNALS
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use QuantityTemplates, only: QuantityTemplate_T
  use VectorsModule, only: ConstructVectorTemplate, VectorTemplate_T
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, SUBTREE

  implicit none
  public
  
  private :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  character (LEN=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module performs the vector template aspects of the construct task

contains ! =====     Public Procedures     =============================

  ! -----------------------------  CreateVecTemplateFromMLSCfInfo  -----
  type (VectorTemplate_T) function CreateVecTemplateFromMLSCfInfo &
    & ( NAME, ROOT, quantityTemplates ) result ( vectorTemplate )

  ! Process the vectorTemplate specification.

    ! Dummy arguments
    integer, intent(in) :: NAME         ! Sub-rosa of name, if any, or zero
    integer, intent(in) :: ROOT         ! Of AST for a vectorTemplate
    type (QuantityTemplate_T), dimension(:) :: quantityTemplates

    ! Local variables
    integer :: I, J           ! Loop inductors
    integer :: nSelections    ! How many selections?
    integer, dimension(:), pointer :: SELECTED
    integer :: SON            ! Son of Root

    ! Executable code

    if ( toggle(gen) ) call &
      & trace_begin ( "ConstructVectorTemplateFromMLSCfInfo", root )

    ! Compute the number of selections
    nSelections = 0
    do i = 2, nsons(root)
      nSelections = nSelections + nsons(subtree(i,root)) - 1
    end do

    call allocate_test ( selected, nSelections, "selected", ModuleName )

    ! Loop through the MLSCF information supplied.  Items are either
    ! lists of quantity template labels, or lists of complete MLS signal
    ! specification strings.

    nSelections = 0
    do i = 2, nsons(root)     ! Skip the "vectorTemplate" name
      son = subtree(i,root)   ! An "assign" vertex of the abstract syntax tree
      select case ( decoration(subtree(1,son)) )
      case ( f_quantities )
        do j = 2, nsons(son)  ! Skip the "quantities" name
          nSelections = nSelections + 1
          ! Get the quantity index that was put into the AST by Construct:
          selected(nSelections) = decoration(decoration(subtree(j,son)))
        end do
      case ( f_signals ) ! ??? Needs work here ???
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          "This version can't handle the SIGNALS field of a vector template" )
      end select
    end do

    ! Now finally construct the vector template

    call ConstructVectorTemplate ( name, quantityTemplates, selected, &
         & vectorTemplate )

    call deallocate_test ( selected, "selected", ModuleName )

    if ( toggle(gen) ) call &
      & trace_end ( "ConstructVectorTemplateFromMLSCfInfo" )

  end function CreateVecTemplateFromMLSCfInfo

!=============================================================================
END MODULE ConstructVectorTemplates
!=============================================================================

!
! $Log$
! Revision 2.1  2000/11/16 02:01:03  vsnyder
! Remove unused variable STATUS.
!
! Revision 2.0  2000/09/05 18:57:02  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

