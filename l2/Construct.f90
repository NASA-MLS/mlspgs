! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Construct                ! The construct module for the MLS L2 sw.
!=============================================================================

  use ConstructQuantityTemplates, only: ConstructMinorFrameQuantity, &
    & CreateQtyTemplateFromMLSCfInfo
  use ConstructVectorTemplates, only: CreateVecTemplateFromMLSCfInfo
  use DUMPER, only: DUMP
  use HGrid, only: AddHGridToDatabase, CreateHGridFromMLSCFInfo, &
    & DestroyHGridDatabase, HGrid_T
  use INIT_TABLES_MODULE, only: S_HGRID, S_QUANTITY, S_TIME, S_VECTORTEMPLATE, &
    & S_VGRID
  use MLSCommon, only: L1BInfo_T, MLSChunk_T, MLSInstrumentNoModules
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error
  use MLSSignalNomenclature
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
  use VGrid, only: AddVGridToDatabase, CreateVGridFromMLSCFInfo, &
    & DestroyVGridDatabase, VGrid_T

  implicit none

  public

  private :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  character(len=130) :: id = &
       "$Id$"
  character (len=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  ! This module performs the `construct' task for the level 2 software.  This
  ! task involves constructing templates for vector quantities, vectors and
  ! matrices.

contains ! =====     Public Procedures     =============================

  ! ---------------------------------------------  MLSL2Construct  -----
  subroutine MLSL2Construct ( root, l1bInfo, chunk, &
       & quantityTemplates, vectorTemplates, mifGeolocation )

  ! This is the `main' subroutine for this module

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Root of the tree for the Construct section
    type (L1BInfo_T), intent(in) :: l1bInfo
    type (MLSChunk_T), intent(in) :: chunk
    type (QuantityTemplate_T), dimension(:), pointer :: quantityTemplates
    type (VectorTemplate_T), dimension(:), pointer :: vectorTemplates
    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation

    ! Local variables
    type (VGrid_T), dimension(:), pointer :: vGrids => NULL()
    type (HGrid_T), dimension(:), pointer :: hGrids => NULL()

    integer :: I                ! Loop counter
    integer :: INSTRUMENTMODULE ! Loop counter
    integer :: KEY              ! S_... from Init_Tables_Module.
    integer :: NAME             ! Sub-rosa index of name
    integer :: SON              ! Son or grandson of Root
    integer :: STATUS           ! Flag
    double precision :: T1, T2  ! for timing
    logical :: TIMING

    ! Executable code
    timing = .false.

    ! First we're going to setup our mifGeolocation quantityTemplates.
    ! These are just two quantity templates containing geolocation
    ! information for the GHz and THz modules.  The software can then
    ! point to these for geolocation information for all minor frame
    ! quantities saving file IO and memory.

    if ( toggle(gen) ) call trace_begin ( "MLSL2Construct", root )

    allocate ( mifGeolocation(MLSInstrumentNoModules), STAT=status )
    if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//"mifGeolocation" )
    do instrumentModule = 1, MLSInstrumentNoModules
      call ConstructMinorFrameQuantity ( l1bInfo, chunk, instrumentModule, &
        & mifGeolocation(instrumentModule) )
    end do

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

      select case( decoration(subtree(1,decoration(subtree(1,key)))) )
      case( s_hgrid )
        call decorate ( key, AddHGridToDatabase ( hGrids, &
          & CreateHGridFromMLSCFInfo ( name, key, l1bInfo, chunk ) ) )
      case ( s_vgrid )
        call decorate ( key, AddVGridToDatabase ( vGrids, &
          & CreateVGridFromMLSCFInfo ( name, key ) ) )
      case ( s_quantity )
        call decorate ( key, AddQuantityTemplateToDatabase ( &
          quantityTemplates, CreateQtyTemplateFromMLSCfInfo ( name, key, &
             & hGrids, vGrids, l1bInfo, chunk ) ) )
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
        call dump ( vgrids )
        call dump ( quantityTemplates )
        call dump ( vectorTemplates )
      end if
      call trace_end ( "MLSL2Construct" )
    end if

    call destroyHGridDatabase ( hGrids )
    call destroyVGridDatabase ( vGrids )
  contains
    subroutine SayTime
      call cpu_time ( t2 )
      call output ( "Timing for MLSL2Construct =" )
      call output ( t2 - t1, advance = 'yes' )
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
    call destroyQuantityTemplateDatabase ( quantityTemplates )
    call destroyQuantityTemplateDatabase ( mifGeolocation )
  end subroutine MLSL2DeConstruct

!=============================================================================
END MODULE Construct
!=============================================================================

!
! $Log$
! Revision 2.3  2000/11/16 02:00:48  vsnyder
! Added timing.
!
! Revision 2.2  2000/09/11 19:23:48  ahanzel
! Removed extra log entries in file.
!
! Revision 2.1  2000/09/08 22:55:56  vsnyder
! Revised to use the tree output by the parser
!

