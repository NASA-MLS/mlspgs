! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module ConstructVectorTemplates ! Construct a template for a vector
!=============================================================================

  use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use DUMP_0, only: DUMP
  use HIGHOUTPUT, only: OUTPUTNAMEDVALUE
  use INIT_TABLES_MODULE, only: F_ADOPT, F_QUANTITIES, &
    & F_REMOVEQUANTITIES, F_REMOVETEMPLATE, F_SOURCE, F_TEMPLATE, &
    & FIELD_FIRST, FIELD_LAST
  use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
  use MLSSETS, only: RELATIVECOMPLEMENT
  use OUTPUT_M, only: OUTPUT
  use QUANTITYTEMPLATES, only: QUANTITYTEMPLATE_T
  use STRING_TABLE, only: GET_STRING
  use TOGGLES, only: GEN, TOGGLE, LEVELS
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE, WHERE
  use VECTORSMODULE, only: CONSTRUCTVECTORTEMPLATE, VECTORTEMPLATE_T, &
    & NULLIFYVECTORTEMPLATE
  use L2PC_M, only: ADOPTVECTORTEMPLATE

  implicit none
  private
  public :: CreateVecTemplateFromMLSCfInfo
  
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  logical, parameter :: DEEBUG = .false.

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
    integer :: I, J                     ! Loop inductors
    integer :: NOQUANTITIES             ! How many selections?
    integer :: NOREMOVEDQUANTITIES      ! How many selections?
    integer :: NOUNIQUEQUANTITIES       ! How many selections are unique?
    integer, dimension(:), pointer :: QUANTITIES
    integer, dimension(:), allocatable :: RELATIVEQUANTITIES
    integer, dimension(:), pointer :: REMOVEDQUANTITIES
    integer :: SON                      ! Son of Root
    integer :: KEY                      ! Son of son
    integer :: ADOPTBIN                 ! Name of l2pc bin to adopt
    integer :: SOURCE                   ! l_rows or l_columns
    logical, dimension(field_first:field_last) :: GOT ! Fields
    character(len=80) :: MESSAGE        ! Possible error message
    logical, parameter :: DIEIFDUP = .false.
        
    ! Executable code

    if ( toggle(gen) .and. levels(gen) > 0 ) call &
      & trace_begin ( "ConstructVectorTemplateFromMLSCfInfo", root )

    call nullifyVectorTemplate ( vectorTemplate ) ! for Sun's still useless compiler
    nullify ( quantities, removedQuantities )
    got = .false.
    noQuantities = 0
    noRemovedQuantities = 0

    ! Loop over the arguments
    ! call outputNamedValue( 'Nsons(root) Construct', nsons(root) )
    do i = 2, nsons(root)     ! Skip the "vectorTemplate" name
      son = subtree(i,root)
      key = subtree(1,son)
      got ( decoration(key) ) = .true.
      select case ( decoration(key) )
      case ( f_adopt )
        adoptBin = sub_rosa(subtree(2,son))
      case ( f_source )
        source = decoration ( subtree(2,son) )
      case ( f_quantities )
        noQuantities = noQuantities + nsons(son) - 1
      case ( f_removeQuantities )
        noRemovedQuantities = noRemovedQuantities + nsons(son) - 1
      case ( f_removeTemplate )
        do j = 2, nsons(son)  ! Skip the "template" name
          call count_quantities ( decoration(subtree(j,son)), noRemovedQuantities )
        end do
      case ( f_template )
        ! call outputNamedValue( 'NoQuantities on entering Construct', noQuantities )
        ! call outputNamedValue( 'Nsons(son) template', nsons(son) )
        do j = 2, nsons(son)  ! Skip the "template" name
          ! call output( 'Entering count_quantities', advance='yes' )
          call count_quantities ( decoration(subtree(j,son)), noQuantities )
        end do
        ! call outputNamedValue( 'NoQuantities on leaving Construct', noQuantities )
      end select
    end do

    if ( got ( f_adopt ) .or. got ( f_source) ) then
      ! Adoption requested
      if ( .not. got ( f_source ) ) call Announce_Error ( key, &
        & 'Must supply source=rows/columns for adoption' )
      vectorTemplate = AdoptVectorTemplate ( adoptBin, name, quantityTemplates, source, message )
      if ( len_trim(message) > 0 ) call Announce_Error ( key, message )
      if ( got ( f_quantities ) ) then
        do j = 1, noQuantities
          if ( all ( quantities(j) /= vectorTemplate%quantities ) ) then
            call get_string ( quantityTemplates(quantities(j))%name, message, strip=.true. )
            call Announce_Error ( key, 'Quantity ' // trim(message) // ' is not in this adopted template' )
          end if
        end do
      end if
      if ( toggle(gen) .and. levels(gen) > 0 ) call &
        & trace_end ( "ConstructVectorTemplateFromMLSCfInfo" )
      return
    else if ( got ( f_quantities ) .or. got ( f_template ) ) then
      ! Not an adoption, just a regular template construction
      call allocate_test ( quantities, noQuantities, "quantities", ModuleName )
      noQuantities = 0
      noUniqueQuantities = 0
      call fill_quantities ( root, quantities )
    else
      ! User didn't give appropriate command
      call Announce_Error ( key, 'Must supply adopt, quantites, source, or template' )
    end if
    ! Do we wish to remove any quantities?
    if ( got ( f_removeQuantities ) .or. got ( f_removeTemplate ) ) then
      if ( DEEBUG ) call outputNamedValue( 'Num of values to be removed', noRemovedQuantities )
      call allocate_test ( removedQuantities, noRemovedQuantities, &
        & "removedQuantities", ModuleName )
      noQuantities = 0
      noUniqueQuantities = 0
      call fill_quantities ( root, removedQuantities, remove=.true. )
      ! Now let's try to remove them
      if ( DEEBUG ) then
        call dump( Quantities, 'Quantities ' )
        call dump( removedQuantities, 'removedQuantities ' )
      endif
      relativeQuantities = RelativeComplement( removedQuantities, Quantities )
      if ( DEEBUG ) call dump( relativeQuantities, 'relativeQuantities ' )
      call ConstructVectorTemplate ( name, &
        & quantityTemplates, relativeQuantities, &
        & vectorTemplate, where=where(root), forWhom=moduleName )
      call deallocate_test ( removedQuantities, "removedQuantities", ModuleName )
      call deallocate_test ( quantities, "quantities", ModuleName )
      call deallocate_test ( relativeQuantities, "relativeQuantities", ModuleName )
    else
      call ConstructVectorTemplate ( name, &
        & quantityTemplates, quantities(1:noUniqueQuantities), &
        & vectorTemplate, where=where(root), forWhom=moduleName )
      call deallocate_test ( quantities, "quantities", ModuleName )
    endif

    if ( toggle(gen) .and. levels(gen) > 0 ) call &
      & trace_end ( "ConstructVectorTemplateFromMLSCfInfo" )

  contains

    recursive subroutine Count_Quantities ( Root, noQuantities )
      integer, intent(in)    :: Root ! of a spec_args in a vectorTemplate
      integer, intent(inout) :: noQuantities ! How many?
      integer :: I, J, Son
      do i = 2, nsons(root)
        son = subtree(i,root)
        select case ( decoration(subtree(1,son)) )
        case ( f_quantities )
          noQuantities = noQuantities + nsons(son) - 1
        case ( f_template )
          ! call outputNamedValue( 'NoQuantities on entering Count', noQuantities )
          do j = 2, nsons(son)  ! Skip the "template" name
            call count_quantities ( decoration(subtree(j,son)), noQuantities )
          end do
          ! call outputNamedValue( 'NoQuantities on leaving Count', noQuantities )
        end select
      end do
    end subroutine Count_Quantities

    recursive subroutine Fill_Quantities ( Root, quantities, remove )
      integer, intent(in)            :: Root ! of a spec_args in a vectorTemplate
      integer, dimension(:), pointer :: quantities
      logical, optional, intent(in)  :: remove
      integer :: I, J, Son
      integer :: qIndex
      logical :: myRemove
      myRemove = .false.
      if ( present(remove) ) myRemove = remove
      if ( myRemove .and. DEEBUG ) call output ( 'Removing quantities', advance='yes' )
      do i = 2, nsons(root)
        son = subtree(i,root)
        if ( myRemove ) then
          if ( any( decoration(subtree(1,son)) == &
            & (/ f_quantities, f_template /) ) ) cycle
        else
          if ( any( decoration(subtree(1,son)) == &
            & (/ f_removequantities, f_removetemplate /) ) ) cycle
        endif
        select case ( decoration(subtree(1,son)) )
        case ( f_quantities, f_removeQuantities )
          do j = 2, nsons(son)
            noQuantities = noQuantities + 1
            ! Get the quantity index that was put into the AST by Construct:
            ! quantities(noQuantities) = decoration(decoration(subtree(j,son)))
            qIndex =  decoration(decoration(subtree(j,son)))
            if ( noUniqueQuantities < 1 ) then
              noUniqueQuantities = noUniqueQuantities + 1
              quantities(noUniqueQuantities) = qIndex
            elseif ( any( qindex == quantities(1:noUniqueQuantities) ) ) then
              if ( DIEIFDUP ) call Announce_Error  ( key, 'Duplicate quantities' )
            else
              noUniqueQuantities = noUniqueQuantities + 1
              quantities(noUniqueQuantities) = qIndex
            endif
          end do ! j = 2, nsons(son)
        case ( f_template, f_removetemplate )
          if ( myRemove .and. DEEBUG ) &
            & call outputnamedValue( 'noQuantities before removeTemplate', &
            & (/ noQuantities, noUniqueQuantities /) )
          do j = 2, nsons(son)  ! Skip the "template" name
          ! The following is a tricky point:
            ! We want to remove the quantities that are part of the template
            ! Note however that when the template was defined, they were
            ! part of a 'quantities" field, not part of a "removeQuantities" field
            call fill_quantities ( decoration(subtree(j,son)), quantities, &
              & remove=.false. )
          end do
          if ( myRemove .and. DEEBUG ) &
            & call outputnamedValue( 'noQuantities after removeTemplate', &
            & (/ noQuantities, noUniqueQuantities /) )
        end select
      end do
    end subroutine Fill_Quantities

  end function CreateVecTemplateFromMLSCfInfo

!=============================================================================

  ! -----------------------------------------------  Announce_Error  -----
  subroutine Announce_Error ( where, message, extra )

    use INTRINSIC, only: LIT_INDICES
    use MORETREE, only: STARTERRORMESSAGE
    use STRING_TABLE, only: DISPLAY_STRING

    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    character (LEN=*), intent(in) :: MESSAGE
    integer, intent(in), optional :: EXTRA

    call startErrorMessage ( where )
    call output ( message )
    if ( present ( extra ) ) call display_string ( lit_indices ( extra ), strip=.true. )
    call output ( '', advance='yes' )
    call MLSMessage ( MLSMSG_Error, ModuleName, 'Problem with vector template construction' )
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

END MODULE ConstructVectorTemplates
!=============================================================================

!
! $Log$
! Revision 2.22  2017/03/10 00:41:31  vsnyder
! Make RELATIVEQUANTITIES allocatable
!
! Revision 2.21  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.20  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.19  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.18  2013/06/12 02:37:14  vsnyder
! Cruft removal
!
! Revision 2.17  2012/06/07 22:43:54  pwagner
! May remove Quantities, vectorTemplates
!
! Revision 2.16  2012/05/24 21:07:38  vsnyder
! Allow any number of template and quantities fields.  Trace out template
! fields recursively until arriving at template fields.  The quantities of
! the created template are the totality of the quantities in quantities fields,
! and in quantities fields of referenced vector templates.
!
! Revision 2.15  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.14  2006/08/05 02:12:27  vsnyder
! Add ForWhom argument to ConstructVectorTemplate
!
! Revision 2.13  2006/07/27 03:52:06  vsnyder
! Pass source_ref to ConstructVectorTemplate
!
! Revision 2.12  2005/06/22 18:57:01  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.11  2004/01/24 01:04:21  livesey
! Added stuff to allow one to adopt quantities
!
! Revision 2.10  2004/01/23 19:08:44  livesey
! Changes / improvements to the adoption.
!
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

