! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module TREE_WALKER

! Traverse the tree output by the parser and checked by the tree checker.
! Perform the actions of the MLS L2 processing in the order indicated.

  use AntennaPatterns_m, only: Destroy_Ant_Patterns_Database
  use Construct, only: MLSL2Construct, MLSL2DeConstruct
  use Dumper, only: Dump
  use Fill, only: MLSL2Fill
  use FilterShapes_m, only: Destroy_Filter_Shapes_Database
  use ForwardModelConfig, only: ForwardModelConfig_T, DestroyFWMConfigDatabase
  use Global_Settings, only: Set_Global_Settings
  use GriddedData, only: GriddedData_T
  use HGrid, only: HGrid_T
  use Init_Tables_Module, only: Z_CHUNKDIVIDE,  Z_CONSTRUCT, Z_FILL, &
    & Z_GLOBALSETTINGS, Z_JOIN, Z_MERGEAPRIORI, Z_MLSSIGNALS, Z_OUTPUT, &
    & Z_READAPRIORI, Z_RETRIEVE, Z_SPECTROSCOPY
  use JOIN, only: MLSL2Join
  use L2AUXData, only: DestroyL2AUXDatabase, L2AUXData_T, Dump_L2AUX
  use L2GPData, only: DestroyL2GPDatabase, L2GPData_T, Dump_L2GP
  use L2ParInfo, only: PARALLEL, CLOSEPARALLEL
  use L2Parallel, only: GETCHUNKFROMMASTER, L2MASTERTASK
  use L2PC_m, only: DestroyL2PCDatabase
  use MatrixModule_1, only: DestroyMatrixDatabase, Matrix_Database_T
  use MLSCommon, only: L1BINFO_T, MLSCHUNK_T, TAI93_RANGE_T
  use MLSSignals_M, only: Bands, DestroyBandDatabase, DestroyModuleDatabase, &
    & DestroyRadiometerDatabase, DestroySignalDatabase, &
    & DestroySpectrometerTypeDatabase, MLSSignals, Modules, Radiometers, &
    & Signals, SpectrometerTypes
  use NCEP_DAO, only: DestroyGridTemplateDatabase
  use Open_Init, only: DestroyL1BInfo, OpenAndInitialize
  use Output_m, only: Output
  use OutputAndClose, only: Output_Close
  use PointingGrid_m, only: Destroy_Pointing_Grid_Database
  use QuantityTemplates, only: QuantityTemplate_T
  use ReadAPriori, only: read_apriori
  use RetrievalModule, only: Retrieve
  use ScanDivide, only: DestroyChunkDatabase, ScanAndDivide
  use SpectroscopyCatalog_m, only: Destroy_Line_Database, &
    & Destroy_SpectCat_Database, Spectroscopy
  use Test_Parse_Signals_m, only: Test_Parse_Signals
  use Toggles, only: GEN, LEVELS, SWITCHES, TOGGLE
  use Trace_m, only: DEPTH, TRACE_BEGIN, TRACE_END
  use Tree, only: DECORATION, NSONS, SUBTREE
  use VectorsModule, only: DestroyVectorDatabase, Vector_T, VectorTemplate_T
  use VGridsDatabase, only: DestroyVGridDatabase, VGrid_T
  use WriteMetadata, only: PCFData_T

  implicit NONE
  private

  public :: WALK_TREE_TO_DO_MLS_L2

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================
  ! -------------------------------------  WALK_TREE_TO_DO_MLS_L2  -----
  subroutine WALK_TREE_TO_DO_MLS_L2 ( ROOT, ERROR_FLAG, FIRST_SECTION )
    integer, intent(in) :: ROOT         ! Root of the abstract syntax tree
    integer, intent(out) :: ERROR_FLAG  ! Nonzero means failure
    integer, intent(in) :: FIRST_SECTION! Index of son of root of first n_cf

    integer :: chunkNo                  ! Index of Chunks
    type (MLSChunk_T), dimension(:), pointer :: Chunks ! of data
    ! Forward model configurations:
    type (ForwardModelConfig_T), dimension(:), &
      & pointer :: ForwardModelConfigDatabase
    type (GriddedData_T), dimension(:), pointer :: GriddedData
    type (HGrid_T), dimension(:), pointer :: HGrids
    integer :: HOWMANY                  ! Nsons(Root)
    integer :: I, J                     ! Loop inductors
    type (L1BInfo_T) :: L1BInfo         ! File handles etc. for L1B dataset
    type (L2AUXData_T), dimension(:), pointer :: L2AUXDatabase
    type (L2GPData_T), dimension(:), pointer  :: L2GPDatabase
    type (PCFData_T) :: L2pcf
    type (Matrix_Database_T), dimension(:), pointer :: Matrices
    type (TAI93_Range_T) :: ProcessingRange  ! Data processing range
    integer :: SON                      ! Son of Root
    type (Vector_T), dimension(:), pointer :: Vectors
    type (VGrid_T), dimension(:), pointer :: VGrids

    ! Arguments for Construct not declared above:
    type (QuantityTemplate_T), dimension(:), pointer :: MifGeolocation
    type (QuantityTemplate_T), dimension(:), pointer :: QtyTemplates
    type (VectorTemplate_T), dimension(:), pointer :: VectorTemplates

    nullify ( chunks, forwardModelConfigDatabase, griddedData, &
      & hGrids, l2auxDatabase, l2gpDatabase, matrices, mifGeolocation, &
      & qtyTemplates, vectors, vectorTemplates, vGrids )

    depth = 0
    if ( toggle(gen) ) call trace_begin ( 'WALK_TREE_TO_DO_MLS_L2', &
      & subtree(first_section,root) )
    call OpenAndInitialize ( processingRange, l1bInfo, l2pcf )

    i = first_section
    howmany = nsons(root)
    do while ( i <= howmany )
      son = subtree(i,root)
      select case ( decoration(subtree(1,son)) ) ! section index
      case ( z_globalsettings )
        call set_global_settings ( son, forwardModelConfigDatabase, vGrids, &
          & l2gpDatabase, l2pcf, processingRange, l1bInfo )
      case ( z_mlsSignals )
        call MLSSignals ( son )
        if ( index(switches,'tps') /= 0 ) call test_parse_signals
      case ( z_spectroscopy )
        call spectroscopy ( son )
      case ( z_readapriori )
      	call read_apriori ( son , l2gpDatabase, l2auxDatabase, griddedData)
      case ( z_mergeapriori )
        ! Merge apriori here
      case ( z_chunkdivide )
        if ( .not. parallel%slave ) then
          call ScanAndDivide ( son, processingRange, l1bInfo, chunks )
        else
          call GetChunkFromMaster ( chunks )
        endif
        if ( toggle(gen) .and. levels(gen) > 0 ) call dump ( chunks )
      case ( z_construct, z_fill, z_join, z_retrieve )
        if ( parallel%master ) then 
          call L2MasterTask ( chunks, l2gpDatabase, l2auxDatabase )
        else
          do chunkNo = lbound(chunks,1), ubound(chunks,1)
            if ( index(switches,'chu') /= 0 ) then
              call output ( " ================ Starting processing for chunk " )
              call output ( chunkNo )
              call output ( " ================ ", advance='yes' )
            endif
            j = i
subtrees:   do while ( j <= howmany )
              son = subtree(j,root)
              select case ( decoration(subtree(1,son)) ) ! section index
              case ( z_construct )
                call MLSL2Construct ( son, l1bInfo, chunks(chunkNo), &
                  & qtyTemplates, vectorTemplates, vGrids, hGrids, &
                  & l2gpDatabase, mifGeolocation )
              case ( z_fill )
                call MLSL2Fill ( son, l1bInfo, griddedData, vectorTemplates, &
                  & vectors, qtyTemplates, matrices, vGrids, l2gpDatabase , &
                  & l2auxDatabase, chunks, chunkNo)
              case ( z_join )
                call MLSL2Join ( son, vectors, l2gpDatabase, &
                  & l2auxDatabase, chunkNo, chunks )
              case ( z_retrieve )
                call retrieve ( son, vectors, matrices, forwardModelConfigDatabase)
              case default
                exit subtrees
              end select
              j = j + 1
            end do subtrees
            ! Now, if we're dealing with more than one chunk destroy stuff
            ! Otherwise, we'll save them as we may need to output them as l2pc files.
            if ( size(chunks) > 1) then
              call MLSL2DeConstruct ( qtyTemplates, vectorTemplates, &
                & mifGeolocation, hGrids )
              call DestroyVectorDatabase ( vectors )
              call DestroyMatrixDatabase ( matrices )
            end if
          end do ! on chunkNo
          i = j - 1 ! one gets added back in at the end of the outer loop
        end if
      case ( z_output ) ! Write out the data
        if ( .not. parallel%slave ) &
          & call Output_Close ( son, l2gpDatabase, l2auxDatabase, matrices, l2pcf,&
          & size(chunks)==1 )

        ! For case where there was one chunk, destroy vectors etc.
        ! This is to guard against destroying stuff needed by l2pc writing
        if ( size(chunks) == 1) then
          call MLSL2DeConstruct ( qtyTemplates, vectorTemplates, &
            & mifGeolocation, hGrids )
          call DestroyVectorDatabase ( vectors )
          call DestroyMatrixDatabase ( matrices )
        end if

        ! Now tidy up any remaining `pointer' data.
        ! processingRange needs no deallocation
        call DestroyL1BInfo ( l1bInfo )
        call DestroyGridTemplateDatabase ( griddedData )
        call DestroyChunkDatabase (chunks )
        if ( index(switches,'l2gp') /= 0 .and. .not. parallel%slave) then
          call Dump_L2GP(l2gpDatabase)
        elseif ( index(switches,'cab') /= 0 .and. .not. parallel%slave) then
          call Dump_L2GP(l2gpDatabase, ColumnsOnly=.true.)
        endif
        call DestroyL2GPDatabase ( l2gpDatabase )
        if ( index(switches,'l2aux') /= 0 .and. .not. parallel%slave) then
          call Dump_L2AUX(l2auxDatabase)
        endif
        call DestroyL2AUXDatabase ( l2auxDatabase )
        ! vectors, vectorTemplates and qtyTemplates destroyed at the
        ! end of each chunk

      end select
      i = i + 1
    end do
    call CloseParallel
    call destroy_ant_patterns_database
    call DestroyL2PCDatabase
    call destroy_filter_shapes_database
    call DestroyFWMConfigDatabase ( forwardModelConfigDatabase )
    call destroy_line_database
    call destroy_pointing_grid_database
    call destroy_spectcat_database
    call DestroyBandDatabase ( Bands )
    call DestroyModuleDatabase ( Modules )
    call DestroyRadiometerDatabase ( Radiometers )
    call DestroySpectrometerTypeDatabase ( SpectrometerTypes )
    call DestroySignalDatabase ( Signals )
    call destroyVGridDatabase ( vGrids )
    error_flag = 0
    if ( toggle(gen) ) call trace_end ( 'WALK_TREE_TO_DO_MLS_L2' )
  end subroutine WALK_TREE_TO_DO_MLS_L2
end module TREE_WALKER

! $Log$
! Revision 2.57  2001/07/11 21:40:21  livesey
! Added -Schu option
!
! Revision 2.56  2001/06/13 20:44:08  livesey
! Moved the CloseParallel higher up, to work around the memory management
! problem (somethine [HDF?] seems to stamp on PointingFrequencyDatabase)
!
! Revision 2.55  2001/06/07 21:58:28  pwagner
! Added Copyright statement
!
! Revision 2.54  2001/05/23 22:00:16  livesey
! Interim version
!
! Revision 2.53  2001/05/23 01:44:35  livesey
! Parallel stuff taking shape
!
! Revision 2.52  2001/05/10 00:43:23  livesey
! Tree walker now owns hGrids
!
! Revision 2.51  2001/05/04 17:12:25  pwagner
! Passes necessary args to global_settings
!
! Revision 2.50  2001/05/03 20:34:08  vsnyder
! Cosmetic changes
!
! Revision 2.49  2001/05/02 23:23:00  livesey
! Added some of the parallel stuff
!
! Revision 2.48  2001/04/28 01:31:46  livesey
! Changes for new l2pc / matrix handling.
!
! Revision 2.47  2001/04/26 20:02:09  livesey
! Made l2pc database a saved array in L2PC_m
!
! Revision 2.46  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.45  2001/04/26 00:08:02  livesey
! Stuff to support reading of l2pc files
!
! Revision 2.44  2001/04/25 21:54:16  livesey
! Added canDoL2PC flag to join
!
! Revision 2.43  2001/04/25 21:52:13  livesey
! Moved DeConstruct to after output for 1 chunk cases.
! This is to protect vectors and matrices stored in l2pcs.
!
! Revision 2.42  2001/04/25 20:34:36  livesey
! Now supports writing of l2pc files
!
! Revision 2.41  2001/04/25 19:31:13  livesey
! Fixed bug, now nullifies l2pcDatabase
!
! Revision 2.40  2001/04/24 23:05:54  vsnyder
! Make 'test_parse_signals' depend on 'switches' containing 'tps'
!
! Revision 2.39  2001/04/24 20:20:02  livesey
! L2PC moved to lib, and renamed
!
! Revision 2.38  2001/04/24 20:05:50  livesey
! New stuff to support joining of l2pc's
!
! Revision 2.37  2001/04/21 01:41:35  vsnyder
! Fix memory leaks
!
! Revision 2.36  2001/04/21 01:26:37  livesey
! Now passes l2gpDatabase to more people
!
! Revision 2.35  2001/04/20 17:12:38  livesey
! Add vGrids argument to fill to support fill from vGrid
!
! Revision 2.34  2001/04/19 23:51:40  pwagner
! Moved anText to become component of PCFData_T
!
! Revision 2.33  2001/04/11 17:47:47  pwagner
! presets anText to null
!
! Revision 2.32  2001/04/10 22:27:47  vsnyder
! Nullify explicitly instead of with <initialization> so as not to give
! pointers the SAVE attribute.  <initialization> is NOT executed on each
! entry to a procedure.
!
! Revision 2.31  2001/04/10 02:46:17  livesey
! Working version, no more FMI/TFMI
!
! Revision 2.30  2001/04/10 00:02:19  vsnyder
! Implement 'matrix' spec in Fill section
!
! Revision 2.29  2001/04/07 01:50:49  vsnyder
! Move some of VGrid to lib/VGridsDatabase.  Move ForwardModelConfig_T and
! some related stuff to fwdmdl/ForwardModelConfig.
!
! Revision 2.28  2001/04/06 20:12:59  vsnyder
! Make 'call test_parse_signal' depend on 'emit' toggle
!
! Revision 2.27  2001/04/04 02:15:12  vsnyder
! Add Spectroscopy section
!
! Revision 2.26  2001/04/03 20:50:45  pwagner
! Added anText to hold PCF file contents
!
! Revision 2.25  2001/04/02 23:41:09  pwagner
! Now keeps l2pcf and transmits as needed
!
! Revision 2.24  2001/03/29 19:13:03  livesey
! Renamed apriorDatabase to griddedData
!
! Revision 2.23  2001/03/28 23:47:48  livesey
! Added arguments to sids etc.
!
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
