module TREE_WALKER

! Traverse the tree output by the parser and checked by the tree checker.
! Perform the actions of the MLS L2 processing in the order indicated.

  use Construct, only: MLSL2Construct, MLSL2DeConstruct
  use DUMPER, only: DUMP
  use FILL, only: MLSL2Fill
  use GLOBAL_SETTINGS, only: SET_GLOBAL_SETTINGS
  use GriddedData, only: DestroyGridTemplateDatabase, GriddedData_T
  use INIT_TABLES_MODULE, only: Z_CHUNKDIVIDE, Z_CONSTRUCT, Z_FILL, &
    & Z_GLOBALSETTINGS, Z_JOIN, Z_MERGEAPRIORI, Z_OUTPUT, Z_READAPRIORI
  use JOIN, only: MLSL2Join
  use L2AUXData, only: DestroyL2AUXDatabase, L2AUXData_T
  use L2GPData, only: DestroyL2GPDatabase, L2GPData_T
  use MLSCommon, only: L1BINFO_T, MLSCHUNK_T, TAI93_RANGE_T
  use ObtainClimatology, only: OBTAIN_CLIM
  use ObtainDAO, only: OBTAIN_DAO
  use ObtainNCEP, only: OBTAIN_NCEP
  use OPEN_INIT, only: DestroyL1BInfo, OpenAndInitialize
  use OutputAndClose, only: Output_Close
  use QuantityTemplates, only: QuantityTemplate_T
  use ScanDivide, only: DestroyChunkDatabase, ScanAndDivide
  use TOGGLES, only: GEN, LEVELS, TOGGLE
  use TRACE_M, only: DEPTH, TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATION, NSONS, SUBTREE
  use TREE_TYPES ! Everything, especially everything beginning with N_
  use VectorsModule, only: DestroyVectorDatabase, Vector_T, VectorTemplate_T

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
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase
    type (L2GPData_T), dimension(:), pointer  :: l2gpDatabase
    type (TAI93_Range_T) :: ProcessingRange  ! Data processing range
    integer :: SON                      ! Son of Root
    type (Vector_T), dimension(:), pointer :: Vectors => NULL()

    ! Arguments for Construct not declared above:
    type (QuantityTemplate_T), dimension(:), pointer :: qtyTemplates => NULL()
    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation => NULL()
    type (VectorTemplate_T), dimension(:), pointer :: vectorTemplates => NULL()

    depth = 0
    if ( toggle(gen) ) call trace_begin ( 'WALK_TREE_TO_DO_MLS_L2', root )
    call OpenAndInitialize ( processingRange, l1bInfo )
    ! For now, the next three are simply done.  Eventually, they should be
    ! triggered by a command
    call Obtain_NCEP ( aprioriData, root )
    call Obtain_DAO ( aprioriData, root )
    call Obtain_Clim ( aprioriData, root )

    i = first_section
    howmany = nsons(root)
    do while ( i <= howmany )
      son = subtree(i,root)
      select case ( decoration(subtree(1,son)) ) ! section index
      case ( z_globalsettings )
        call set_global_settings ( son )
      case ( z_readapriori )
        ! Read apriori here
      case ( z_mergeapriori )
        ! Merge apriori here
      case ( z_chunkdivide )
        call ScanAndDivide ( son, processingRange, l1bInfo, chunks )
        if ( toggle(gen) .and. levels(gen) > 0 ) call dump ( chunks )
      case ( z_construct, z_fill, z_join )
        do chunkNo = 1, size(chunks)
          j = i
subtrees: do while ( j <= howmany )
            son = subtree(j,root)
            select case ( decoration(subtree(1,son)) ) ! section index
            case ( z_construct )
              call MLSL2Construct ( son, l1bInfo, chunks(chunkNo), &
                & qtyTemplates, vectorTemplates, mifGeolocation )
            case ( z_fill )
              call MLSL2Fill ( son, l1bInfo, aprioriData, vectorTemplates, &
                & vectors, qtyTemplates, l2gpDatabase )
            case ( z_join )
              call MLSL2Join ( son, vectors, l2gpDatabase, l2auxDatabase, &
                & qtyTemplates, chunks, chunkNo )
!           case ( z_retrieve )
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
        ! vectors, vectorTemplates and qtyTemplates destroyed at the
        ! end of each chunk
        call DestroyL2GPDatabase ( l2gpDatabase )
        call DestroyL2AUXDatabase ( l2auxDatabase )

      end select
      i = i + 1
    end do
    error_flag = 0
    if ( toggle(gen) ) call trace_end ( 'WALK_TREE_TO_DO_MLS_L2' )
  end subroutine WALK_TREE_TO_DO_MLS_L2
end module TREE_WALKER

! $Log$
