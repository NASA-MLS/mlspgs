! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Construct                ! The construct module for the MLS L2 sw.
!=============================================================================

  use Allocate_Deallocate, only: Deallocate_test
  use ConstructQuantityTemplates, only: ConstructMinorFrameQuantity, &
    & CreateQtyTemplateFromMLSCfInfo, ForgeMinorFrames
  use ConstructVectorTemplates, only: CreateVecTemplateFromMLSCfInfo
  use Dumper, only: Dump
  use HGrid, only: AddHGridToDatabase, CreateHGridFromMLSCFInfo, &
    & DestroyHGridDatabase, HGrid_T
  use INIT_TABLES_MODULE, only: S_FORGE, S_HGRID, S_QUANTITY, S_TIME, S_VECTORTEMPLATE
  use MLSCommon, only: L1BInfo_T, MLSChunk_T
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error
  use MLSSignals_m, only: Modules
  use MoreTree, only: Get_Spec_ID
  use OUTPUT_M, only: OUTPUT
  use QuantityTemplates, only: AddQuantityTemplateToDatabase, &
    & DestroyQuantityTemplateDatabase, QuantityTemplate_T
  use TOGGLES, only: GEN, LEVELS, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, NULL_TREE, SUB_ROSA, &
    & SUBTREE
  use TREE_TYPES, only: N_NAMED
  use VectorsModule, only: AddVectorTemplateToDatabase, &
    & DestroyVectorTemplateDatabase, Dump, VectorTemplate_T
  use VGridsDatabase, only: VGrid_T
  use Intrinsic, ONLY: L_None
  use String_Table, ONLY: GET_STRING
  use Init_tables_module, ONLY: LIT_INDICES
  use L2GPData, only: L2GPDATA_T

  implicit none

  public

  private :: Id, IdParm, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = &
    & "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  ! This module performs the `construct' task for the level 2 software.  This
  ! task involves constructing templates for vector quantities, vectors and
  ! matrices.

contains ! =====     Public Procedures     =============================

  ! ---------------------------------------------  MLSL2Construct  -----
  subroutine MLSL2Construct ( root, l1bInfo, chunk, &
       & quantityTemplates, vectorTemplates, VGrids, l2gpDatabase, mifGeolocation )

  ! This is the `main' subroutine for this module

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Root of the tree for the Construct section
    type (L1BInfo_T), intent(in) :: l1bInfo
    type (MLSChunk_T), intent(in) :: chunk
    type (QuantityTemplate_T), dimension(:), pointer :: quantityTemplates
    type (VectorTemplate_T), dimension(:), pointer :: vectorTemplates
    type (VGrid_T), dimension(:), pointer :: vGrids
    type (L2GPData_T), dimension(:), pointer :: L2GPDatabase
    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation

    ! Local variables
    type (HGrid_T), dimension(:), pointer :: HGrids

    integer :: I                ! Loop counter
    integer :: InstrumentModuleIndex ! Loop counter
    integer :: KEY              ! S_... from Init_Tables_Module.
    integer :: NAME             ! Sub-rosa index of name
    integer :: SON              ! Son or grandson of Root
    integer :: STATUS           ! Flag
    REAL :: T1, T2              ! for timing
    logical :: TIMING

    ! Executable code
    timing = .false.

    ! First we're going to setup our mifGeolocation quantityTemplates.
    ! These are just two quantity templates containing geolocation
    ! information for the GHz and THz modules.  The software can then
    ! point to these for geolocation information for all minor frame
    ! quantities saving file IO and memory.

    if ( toggle(gen) ) call trace_begin ( "MLSL2Construct", root )

    nullify ( hGrids )

    allocate ( mifGeolocation(size(modules)), STAT=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//"mifGeolocation" )
    
    ! Now try to fill it if we have any L1BFiles
    if (l1bInfo%l1boaID /= 0 ) then
      do instrumentModuleIndex = 1, size(modules)
        call ConstructMinorFrameQuantity ( l1bInfo, chunk, &
          & instrumentModuleIndex, mifGeolocation(instrumentModuleIndex) )
      end do
    else
      mifGeolocation%noSurfs = 0
      mifGeolocation%noInstances = 0
    endif

    ! The rest is fairly simple really.  We just loop over the mlscf 
    ! instructions and hand them off to people

    do i = 2, nsons(root)-1 ! Skip the section name at begin and end
      son = subtree(i,root)
      if ( node_id(son) == n_named ) then ! Is spec labeled?
        key = subtree(2,son)
        name = sub_rosa(subtree(1,son))
      else ! Son is n_spec_args
        key = son
        name = 0
      end if

      ! Node_id(key) is now n_spec_args.
      
      select case( get_spec_id(key) )
      case ( s_forge )
        call ForgeMinorFrames ( key, mifGeolocation )
      case( s_hgrid )
        call decorate ( key, AddHGridToDatabase ( hGrids, &
          & CreateHGridFromMLSCFInfo ( name, key, l1bInfo, l2gpDatabase, chunk ) ) )
      case ( s_quantity )
        call decorate ( key, AddQuantityTemplateToDatabase ( &
          & quantityTemplates, CreateQtyTemplateFromMLSCfInfo ( name, key, &
            & hGrids, vGrids, l1bInfo, chunk, mifGeolocation ) ) )
      case ( s_vectortemplate )
        call decorate ( key, AddVectorTemplateToDatabase ( vectorTemplates, &
          & CreateVecTemplateFromMLSCfInfo ( name, key, quantityTemplates ) ) )
      case ( s_time )
        if ( timing ) then
          call sayTime
        else
          call cpu_time ( t1 )
          timing = .true.
        end if
      case default ! Can't get here if tree_checker worked correctly
      end select
    end do

    if ( toggle(gen) ) then
      if (  levels(gen) > 0 ) then
        call dump ( hgrids )
        call dump ( quantityTemplates, details=levels(gen)-1 )
        call dump ( vectorTemplates, details=levels(gen)-1 )
      end if
      call trace_end ( "MLSL2Construct" )
    end if

    call destroyHGridDatabase ( hGrids )
    if ( timing ) call sayTIme

  contains
    subroutine SayTime
      call cpu_time ( t2 )
      call output ( "Timing for MLSL2Construct =" )
      call output ( DBLE(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime
  end subroutine MLSL2Construct

  ! -------------------------------------------  MLSL2DeConstruct  -----
  subroutine MLSL2DeConstruct ( quantityTemplates, vectorTemplates, &
    &                           mifGeolocation )

  ! DeConstruct the Quantity and Vector template databases.

    type (QuantityTemplate_T), dimension(:), pointer :: quantityTemplates
    type (VectorTemplate_T), dimension(:), pointer :: vectorTemplates
    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation

    call destroyVectorTemplateDatabase ( vectorTemplates )
    call destroyQuantityTemplateDatabase ( quantityTemplates, ignoreMinorFrame=.true. )
    call destroyQuantityTemplateDatabase ( mifGeolocation )
  end subroutine MLSL2DeConstruct

!=============================================================================
END MODULE Construct
!=============================================================================

!
! $Log$
! Revision 2.19  2001/04/23 23:53:03  livesey
! Tidied up deallocation of minor frame quantities
!
! Revision 2.18  2001/04/23 23:24:55  livesey
! Changed l2gpDatabase to pointer
!
! Revision 2.17  2001/04/21 01:25:23  livesey
! New version, can construct h/v grids from l2gp
!
! Revision 2.16  2001/04/20 23:11:26  livesey
! Added `Forge' stuff
!
! Revision 2.15  2001/04/10 23:44:44  vsnyder
! Improve 'dump'
!
! Revision 2.14  2001/04/10 22:27:47  vsnyder
! Nullify explicitly instead of with <initialization> so as not to give
! pointers the SAVE attribute.  <initialization> is NOT executed on each
! entry to a procedure.
!
! Revision 2.13  2001/04/07 01:50:48  vsnyder
! Move some of VGrid to lib/VGridsDatabase.  Move ForwardModelConfig_T and
! some related stuff to fwdmdl/ForwardModelConfig.
!
! Revision 2.12  2001/03/28 18:09:24  vsnyder
! Cosmetic changes
!
! Revision 2.11  2001/03/28 01:24:55  vsnyder
! Move vGrid from construct section to global settings section
!
! Revision 2.10  2001/03/15 21:09:52  vsnyder
! Use 'get_spec_id' from More_Tree
!
! Revision 2.9  2001/03/15 21:03:46  vsnyder
! Cross-references between databases are by database index, not tree index
!
! Revision 2.8  2001/03/03 00:07:51  livesey
! New mifGeolocation stuff
!
! Revision 2.7  2001/03/02 01:28:12  livesey
! Uses new MLSSignals
!
! Revision 2.6  2001/02/09 00:38:22  livesey
! Various updates
!
! Revision 2.5  2001/01/03 17:51:05  pwagner
! Changed types of t1, t2 to real
!
! Revision 2.4  2000/11/16 02:15:14  vsnyder
! More work on timing.
!
! Revision 2.3  2000/11/16 02:00:48  vsnyder
! Added timing.
!
! Revision 2.2  2000/09/11 19:23:48  ahanzel
! Removed extra log entries in file.
!
! Revision 2.1  2000/09/08 22:55:56  vsnyder
! Revised to use the tree output by the parser
!

