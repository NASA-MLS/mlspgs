! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ConstructVectorTemplates ! Construct a template for a vector
!=============================================================================

  use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
  use INIT_TABLES_MODULE, only: F_QUANTITIES
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use Output_M, only: Output
  use QuantityTemplates, only: QuantityTemplate_T
  use String_Table, only: Display_String
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, SOURCE_REF, SUB_ROSA, &
    & SUBTREE
  use VectorsModule, only: ConstructVectorTemplate, VectorTemplate_T

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
    integer :: I, J, K        ! Loop inductors
    integer :: nSelections    ! How many selections?
    integer, dimension(:), pointer :: SELECTED
    integer :: SON            ! Son of Root
    integer :: SOURCE         ! 256*line + column of erroneous input

    ! Executable code

    if ( toggle(gen) ) call &
      & trace_begin ( "ConstructVectorTemplateFromMLSCfInfo", root )

    nullify ( selected )

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
      ! Currently only one field, but we'll leave the case statement in to be sure
      select case ( decoration(subtree(1,son)) )
      case ( f_quantities )
        do j = 2, nsons(son)  ! Skip the "quantities" name
          nSelections = nSelections + 1
          ! Get the quantity index that was put into the AST by Construct:
          selected(nSelections) = decoration(decoration(subtree(j,son)))
          ! Check for duplicate quantities
          do k = 1, nSelections - 1
            if ( selected(k) == selected(nSelections) ) then
              source = source_ref( subtree(j,son) )
              call output ( 'At line '  )
              call output ( mod(source,256) )
              call output ( ', column ' )
              call output ( source/256 )
              call output ( ', the quantity ' )
              call display_string( sub_rosa(subtree(j,son)) )
              call output ( ' duplicates a previous one.', advance='yes' )
              call MLSMessage ( MLSMSG_Error, ModuleName, &
                "Duplicate quantity specified in vector template" )
            end if
          end do ! k = 1, nSelections - 1
        end do ! j = 2, nsons(son)
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
! Revision 2.5  2001/10/15 22:05:20  livesey
! Got rid of the signals stuff which was never implemented
!
! Revision 2.4  2001/04/10 22:27:47  vsnyder
! Nullify explicitly instead of with <initialization> so as not to give
! pointers the SAVE attribute.  <initialization> is NOT executed on each
! entry to a procedure.
!
! Revision 2.3  2001/02/22 21:58:42  livesey
! Nullified a pointer
!
! Revision 2.2  2000/12/19 20:14:57  vsnyder
! Add test for duplicate quantities.
!
! Revision 2.1  2000/11/16 02:01:03  vsnyder
! Remove unused variable STATUS.
!
! Revision 2.0  2000/09/05 18:57:02  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!

