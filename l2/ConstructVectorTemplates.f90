! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ConstructVectorTemplates ! Construct a template for a vector
!=============================================================================

  use Allocate_Deallocate, only: Allocate_Test, DeAllocate_Test
  use INIT_TABLES_MODULE, only: F_QUANTITIES, FIELD_FIRST, FIELD_LAST, &
    & F_ADOPTROWS, F_ADOPTCOLUMNS, L_ROWS, L_COLUMNS
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use Output_M, only: Output
  use QuantityTemplates, only: QuantityTemplate_T
  use String_Table, only: Display_String
  use TOGGLES, only: GEN, TOGGLE, LEVELS
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, SOURCE_REF, SUB_ROSA, &
    & SUBTREE
  use VectorsModule, only: ConstructVectorTemplate, VectorTemplate_T, &
    & NullifyVectorTemplate
  use L2PC_m, only: ADOPTVECTORTEMPLATE

  implicit none
  public
  
  private :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  character (LEN=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here 
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
    type (QuantityTemplate_T), dimension(:), pointer :: quantityTemplates

    ! Local variables
    integer :: I, J, K        ! Loop inductors
    integer :: NOQUANTITIES    ! How many selections?
    integer, dimension(:), pointer :: QUANTITIES
    integer :: SON            ! Son of Root
    integer :: KEY            ! Son of son
    logical, dimension(field_first:field_last) :: GOT ! Fields
        
    ! Executable code

    if ( toggle(gen) .and. levels(gen) > 0 ) call &
      & trace_begin ( "ConstructVectorTemplateFromMLSCfInfo", root )

    call nullifyVectorTemplate ( vectorTemplate ) ! for Sun's still useless compiler
    nullify ( quantities )
    got = .false.
    noQuantities = 0

    ! Loop over the arguments
    do i = 2, nsons(root)     ! Skip the "vectorTemplate" name
      son = subtree(i,root)
      key = subtree(1,son)
      got ( decoration(key) ) = .true.
      select case ( decoration(key) )
      case ( f_adoptColumns )
        vectorTemplate = AdoptVectorTemplate ( sub_rosa(subtree(2,son)), quantityTemplates, &
          & source=l_columns )
      case ( f_adoptRows )
        vectorTemplate = AdoptVectorTemplate ( sub_rosa(subtree(2,son)), quantityTemplates, &
          & source=l_rows )
      case ( f_quantities )
        noQuantities = nsons(son) - 1
        call allocate_test ( quantities, noQuantities, "quantities", ModuleName )
        do j = 2, nsons(son)  ! Skip the "quantities" name
          ! Get the quantity index that was put into the AST by Construct:
          quantities(j-1) = decoration(decoration(subtree(j,son)))
          ! Check for duplicate quantities
          do k = 1, j-2
            if ( quantities(k) == quantities(j-1) ) &
              & call Announce_Error  ( key, 'Duplicate quantities' )
          end do ! k = 1, noQuantities - 1
        end do ! j = 2, nsons(son)
      end select
    end do

    if ( (  got ( f_adoptColumns ) .or. got ( f_adoptRows ) ) .and. &
      & .not. associated ( vectorTemplate%quantities ) ) call Announce_Error ( key, &
      & 'No such l2pc bin to adopt' )

    if ( got ( f_quantities ) ) then
      if ( got ( f_adoptColumns ) .or. got ( f_adoptRows ) ) &
        & call Announce_Error ( key, 'Cannot supply both quantities and adoptRows/Columns' )
      call ConstructVectorTemplate ( name, quantityTemplates, quantities, &
        & vectorTemplate )
      call deallocate_test ( quantities, "quantities", ModuleName )
    end if

    if ( toggle(gen) .and. levels(gen) > 0 ) call &
      & trace_end ( "ConstructVectorTemplateFromMLSCfInfo" )

    if ( got ( f_adoptColumns ) .and. got ( f_adoptRows ) ) &
      & call Announce_Error ( key, 'Cannot supply both adoptColumns and adoptRows' )

  end function CreateVecTemplateFromMLSCfInfo

!=============================================================================

  ! -----------------------------------------------  Announce_Error  -----
  subroutine Announce_Error ( where, message, extra )

    use LEXER_CORE, only: PRINT_SOURCE
    use OUTPUT_M, only: OUTPUT
    use TREE, only: SOURCE_REF
    use Intrinsic, only: LIT_INDICES
    use String_Table, only: DISPLAY_STRING
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR

    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    character (LEN=*), intent(in) :: MESSAGE
    integer, intent(in), optional :: EXTRA

    call output ( '***** At ' )
    if ( where > 0 ) then
      call print_source ( source_ref(where) )
    else
      call output ( '(no lcf tree available)' )
    end if
    call output ( ': ' )
    call output ( message )
    if ( present ( extra ) ) call display_string ( lit_indices ( extra ), strip=.true. )
    call output ( '', advance='yes' )
    call MLSMessage ( MLSMSG_Error, ModuleName, 'Problem with vector template construction' )
  end subroutine Announce_Error

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

END MODULE ConstructVectorTemplates
!=============================================================================

!
! $Log$
! Revision 2.9  2004/01/23 05:47:38  livesey
! Added the adoption stuff
!
! Revision 2.8  2002/11/22 12:16:44  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.7  2002/10/08 17:36:20  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.6  2002/09/25 20:08:05  livesey
! Made -g less verbose
!
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

