module TREE_WALKER

! Traverse the tree output by the parser and checked by the tree checker.
! Perform the actions of the MLS L2 processing in the order indicated.

! use Test_Parse_Signals_m, only: Test_Parse_Signals ! Uncomment to test Parse_Signals
  use Construct, only: MLSL2Construct, MLSL2DeConstruct
  use DUMPER, only: DUMP
  use FILL, only: MLSL2Fill
  use ForwardModelInterface, only: ForwardModelConfig_T, DestroyFWMConfigDatabase
  use GLOBAL_SETTINGS, only: SET_GLOBAL_SETTINGS
  use GriddedData, only: GriddedData_T
  use INIT_TABLES_MODULE, only: Field_Indices, Spec_Indices, Z_CHUNKDIVIDE, &
    & Z_CONSTRUCT, Z_FILL, Z_GLOBALSETTINGS, Z_JOIN, Z_MERGEAPRIORI,& 
    & Z_MLSSIGNALS, Z_OUTPUT, Z_READAPRIORI, Z_RETRIEVE
  use JOIN, only: MLSL2Join
  !??? The next USE statement is Temporary for l2load:
  use L2_TEST_STRUCTURES_M, only: FWD_MDL_CONFIG, FWD_MDL_INFO, &
    & TEMPORARY_FWD_MDL_INFO
  use L2AUXData, only: DestroyL2AUXDatabase, L2AUXData_T
  use L2GPData, only: DestroyL2GPDatabase, L2GPData_T
  use MatrixModule_1, only: Matrix_Database_T
  use MLSCommon, only: L1BINFO_T, MLSCHUNK_T, TAI93_RANGE_T
  use MLSSignals_M, only: MLSSignals
!  use OPEN_INIT, only: DestroyL1BInfo, OpenAndInitialize, read_apriori
  use ncep_dao, only: DestroyGridTemplateDatabase
  use OPEN_INIT, only: DestroyL1BInfo, OpenAndInitialize
  use ReadAPriori, only: read_apriori
  use OutputAndClose, only: Output_Close
  use QuantityTemplates, only: QuantityTemplate_T
  use RetrievalModule, only: Retrieve
  use ScanDivide, only: DestroyChunkDatabase, ScanAndDivide
  use TOGGLES, only: GEN, LEVELS, TOGGLE
  use TRACE_M, only: DEPTH, TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATION, NSONS, SUBTREE
  use TREE_TYPES ! Everything, especially everything beginning with N_
  use VectorsModule, only: DestroyVectorDatabase, Vector_T, VectorTemplate_T
  use VGrid, only: DestroyVGridDatabase, VGrid_T

  implicit NONE
  private

  public :: WALK_TREE_TO_DO_MLS_L2

!---------------------------- RCS Ident Info ---------------------------
  character (len=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!-----------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================
  ! -------------------------------------  WALK_TREE_TO_DO_MLS_L2  -----
  subroutine WALK_TREE_TO_DO_MLS_L2 ( ROOT, ERROR_FLAG, FIRST_SECTION )
    integer, intent(in) :: ROOT         ! Root of the abstract syntax tree
    integer, intent(out) :: ERROR_FLAG  ! Nonzero means failure
    integer, intent(in) :: FIRST_SECTION! Index of son of root of first n_cf

    type (GriddedData_T), dimension(:), pointer :: aprioriData => NULL() 
    integer :: chunkNo                  ! Index of Chunks
    type (MLSChunk_T), dimension(:), pointer :: CHUNKS => NULL() ! of data
    integer :: HOWMANY                  ! Nsons(Root)
    integer :: I, J                     ! Loop inductors
    type (L1BInfo_T) :: L1BInfo         ! File handles etc. for L1B dataset
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase => NULL()
    type (L2GPData_T), dimension(:), pointer  :: l2gpDatabase => NULL()
    type (Matrix_Database_T), dimension(:), pointer :: Matrices => NULL()
    type (TAI93_Range_T) :: ProcessingRange  ! Data processing range
    integer :: SON                      ! Son of Root
    type (Vector_T), dimension(:), pointer :: Vectors => NULL()
    type (VGrid_T), dimension(:), pointer :: VGrids => NULL()

    ! Forward model configurations
    type (ForwardModelConfig_T), dimension(:), &
      & pointer :: ForwardModelConfigDatabase => NULL()

    ! Arguments for Construct not declared above:
    type (QuantityTemplate_T), dimension(:), pointer :: qtyTemplates => NULL()
    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation => NULL()
    type (VectorTemplate_T), dimension(:), pointer :: vectorTemplates => NULL()

!??? Begin temporary stuff to start up the forward model
  type(fwd_mdl_config) :: FMC
  type(fwd_mdl_info), dimension(:), pointer :: FMI => NULL()
  type(temporary_fwd_mdl_info), dimension(:), pointer :: TFMI => NULL()
!??? End of temporary stuff to start up the forward model

    depth = 0
    if ( toggle(gen) ) call trace_begin ( 'WALK_TREE_TO_DO_MLS_L2', &
      & subtree(first_section,root) )
    call OpenAndInitialize ( processingRange, l1bInfo )

    i = first_section
    howmany = nsons(root)
    do while ( i <= howmany )
      son = subtree(i,root)
      select case ( decoration(subtree(1,son)) ) ! section index
      case ( z_globalsettings )
!       call set_global_settings ( son, forwardModelConfigDatabase, vGrids ) !??? Restore when l2load isn't needed
        call set_global_settings ( son, forwardModelConfigDatabase, vGrids, &
          & fmc, fmi, tfmi ) !??? This line is temporary for l2load
      case ( z_mlsSignals )
        call MLSSignals ( son, field_indices )
!       call test_parse_signals  ! Uncomment this to test Parse_Signals
      case ( z_readapriori )
        ! Read apriori here
      	CALL read_apriori ( son , l2gpDatabase, l2auxDatabase, aprioriData)
      case ( z_mergeapriori )
        ! Merge apriori here
      case ( z_chunkdivide )
        call ScanAndDivide ( son, processingRange, l1bInfo, chunks )
        if ( toggle(gen) .and. levels(gen) > 0 ) call dump ( chunks )
      case ( z_construct, z_fill, z_join, z_retrieve )
        do chunkNo = 1, size(chunks)
          j = i
subtrees: do while ( j <= howmany )
            son = subtree(j,root)
            select case ( decoration(subtree(1,son)) ) ! section index
            case ( z_construct )
              call MLSL2Construct ( son, l1bInfo, chunks(chunkNo), &
                & qtyTemplates, vectorTemplates, vGrids, mifGeolocation )
            case ( z_fill )
              call MLSL2Fill ( son, l1bInfo, aprioriData, vectorTemplates, &
                & vectors, qtyTemplates, l2gpDatabase , l2auxDatabase, &
                & chunks, chunkNo)
            case ( z_join )
              call MLSL2Join ( son, vectors, l2gpDatabase, l2auxDatabase, chunkNo, chunks )
            case ( z_retrieve )
!             call retrieve ( son, vectors, matrices, forwardModelInfo ) !??? Restore when l2load isn't needed
              call retrieve ( son, vectors, matrices, forwardModelConfigDatabase, &
                & fmc, fmi, tfmi ) !??? This line is temporary for l2load
            case default
          exit subtrees
            end select
            j = j + 1
          end do subtrees
          call MLSL2DeConstruct ( qtyTemplates, vectorTemplates, &
            & mifGeolocation )
          call DestroyVectorDatabase ( vectors )
        end do ! on chunkNo
        i = j - 1 ! one gets added back in at the end of the outer loop
      case ( z_output ) ! Write out the data
        call Output_Close ( son, l2gpDatabase, l2auxDatabase )

        ! Now tidy up any remaining `pointer' data.
        ! processingRange needs no deallocation
        call DestroyL1BInfo ( l1bInfo )
        call DestroyGridTemplateDatabase ( aprioriData )
        call DestroyChunkDatabase (chunks )
        call DestroyL2GPDatabase ( l2gpDatabase )
        call DestroyL2AUXDatabase ( l2auxDatabase )
        call DestroyFWMConfigDatabase ( forwardModelConfigDatabase )
        ! vectors, vectorTemplates and qtyTemplates destroyed at the
        ! end of each chunk

      end select
      i = i + 1
    end do
    call destroyVGridDatabase ( vGrids )
    error_flag = 0
    if ( toggle(gen) ) call trace_end ( 'WALK_TREE_TO_DO_MLS_L2' )
  end subroutine WALK_TREE_TO_DO_MLS_L2
end module TREE_WALKER

! $Log$
! Revision 2.22  2001/03/28 01:24:55  vsnyder
! Move vGrid from construct section to global settings section
!
! Revision 2.21  2001/03/17 03:30:25  vsnyder
! Remove FMI and TFMI from the call to set_global_settings, since it no
! longer uses them.
!
! Revision 2.20  2001/03/17 00:45:53  livesey
! Added forwardModelConfigDatabase
!
! Revision 2.19  2001/03/15 23:26:09  livesey
! Added chunks to call to MLSL2Join
!
! Revision 2.18  2001/03/15 21:20:54  pwagner
! Split between GriddedData and ncep_dao modules
!
! Revision 2.17  2001/03/15 21:05:23  vsnyder
! Set up to test Parse_Signal
!
! Revision 2.16  2001/03/14 02:04:53  vsnyder
! Moved MLSSignals_m to mlspgs/lib
!
! Revision 2.15  2001/03/09 19:57:48  vsnyder
! Add 'z_retrieve' to a case branch
!
! Revision 2.14  2001/03/08 18:21:11  vsnyder
! Even more stuff for L2_Load
!
! Revision 2.13  2001/03/08 03:23:10  vsnyder
! More stuff to work with L2_Load
!
! Revision 2.12  2001/03/07 22:46:05  vsnyder
! Add temporary stuff for Zvi's "l2_load", which will wither away.
!
! Revision 2.11  2001/03/03 00:13:30  pwagner
! Gets read_apriori from ReadAPriori module
!
! Revision 2.10  2001/02/28 01:17:57  livesey
! Removed obtain_ncep etc. These will later be in lower down modules
!
! Revision 2.9  2001/02/27 17:38:07  livesey
! Tidied up arguments to MLSL2Join
!
! Revision 2.8  2001/02/23 02:51:44  vsnyder
! Improve progress messages triggered by -g option
!
! Revision 2.7  2001/02/08 01:40:53  vsnyder
! Don't know what I have done, but "cvs update" said M instead of U
!
! Revision 2.6  2001/01/03 17:48:43  pwagner
! added chunk args to call to Fill
!
! Revision 2.5  2000/12/05 00:41:50  pwagner
! Added L2AUXDatabase arg in call to MLSL2Fill
!
! Revision 2.4  2000/11/30 00:21:08  pwagner
! passes l2*databses to read a priori
!
! Revision 2.3  2000/11/29 00:27:54  pwagner
! Began changes to open old l2gp
!
! Revision 2.2  2000/10/10 00:37:46  vsnyder
! Added $Log for CVS at end.
!
