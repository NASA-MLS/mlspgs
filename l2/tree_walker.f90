! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module TREE_WALKER

! Traverse the tree output by the parser and checked by the tree checker.
! Perform the actions of the MLS L2 processing in the order indicated.


  implicit NONE
  private

  public :: WALK_TREE_TO_DO_MLS_L2

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================
  ! -------------------------------------  WALK_TREE_TO_DO_MLS_L2  -----
  subroutine WALK_TREE_TO_DO_MLS_L2 ( ROOT, ERROR_FLAG, FIRST_SECTION, &
    & COUNTCHUNKS, SINGLECHUNK, LASTCHUNKIN, FILEDATABASE )

    use AntennaPatterns_m, only: Destroy_Ant_Patterns_Database
    use ChunkDivide_m, only: ChunkDivide, DestroyChunkDatabase, &
      & ReduceChunkDatabase
    use Construct, only: MLSL2Construct, MLSL2DeConstruct, &
      & ConstructMIFGeolocation
    use DirectWrite_m, only: DirectData_T, DestroyDirectDatabase
    use Dumper, only: Dump
    use EmpiricalGeometry, only: ForgetOptimumLon0
    use FGrid, only: FGrid_T, DestroyFGridDatabase
    use Fill, only: MLSL2Fill
    use FilterShapes_m, only: Destroy_Filter_Shapes_Database
    use ForwardModelConfig, only: ForwardModelConfig_T, &
      & DestroyFWMConfigDatabase, &
      & StripForwardModelConfigDatabase
    use ForwardModelSupport, only: printForwardModelTiming, &
      & resetForwardModelTiming
    use Global_Settings, only: Set_Global_Settings
    use GriddedData, only: GriddedData_T, DestroyGriddedDataDatabase, Dump
    use HGridsDatabase, only: HGrid_T
    use HGrid, only: COMPUTEALLHGRIDOFFSETS
    use Init_Tables_Module, only: L_CHISQCHAN, L_CHISQMMAF, L_CHISQMMIF, &
      & Z_CHUNKDIVIDE,  Z_CONSTRUCT, Z_FILL, Z_GLOBALSETTINGS, Z_JOIN, &
      & Z_MERGEGRIDS, Z_MLSSIGNALS, Z_OUTPUT, Z_READAPRIORI, Z_RETRIEVE, &
      & Z_SPECTROSCOPY
    use JOIN, only: MLSL2Join
    use L2AUXData, only: DestroyL2AUXDatabase, L2AUXData_T, Dump
    use L2FWMParallel, only: L2FWMSlaveTask, LaunchFWMSlaves
    use L2GPData, only: DestroyL2GPDatabase, L2GPData_T, Dump
    use L2Parallel, only: GETCHUNKINFOFROMMASTER, L2MASTERTASK
    use L2ParInfo, only: PARALLEL, CLOSEPARALLEL
    use L2PC_m, only: DestroyL2PCDatabase, DestroyBinSelectorDatabase
    use LinearizedForwardModel_m, only: FLUSHLOCKEDBINS
    use MACHINE, only: MLS_GC_NOW, MLS_HOWMANY_GC
    use MatrixModule_1, only: DestroyMatrixDatabase, Matrix_Database_T
    use MergeGridsModule, only: MergeGrids
    use MLSCommon, only: L1BINFO_T, MLSCHUNK_T, TAI93_RANGE_T, MLSFile_T
    ! use MLSFiles, only: MLSFile_T
    use MLSL2Options, only: CHECKPATHS
    use MLSMessageModule, only: MLSMessage, MLSMSG_Info, MLSMSG_Error
    use MLSSignals_M, only: Bands, DestroyBandDatabase, DestroyModuleDatabase, &
      & DestroyRadiometerDatabase, DestroySignalDatabase, &
      & DestroySpectrometerTypeDatabase, MLSSignals, Modules, Radiometers, &
      & Signals, SpectrometerTypes
    use MLSL2Timings, only: add_to_section_timing, TOTAL_TIMES
    use Open_Init, only: DestroyL1BInfo, OpenAndInitialize
    use OutputAndClose, only: Output_Close
    use Output_m, only: BLANKS, Output
    use PointingGrid_m, only: Destroy_Pointing_Grid_Database
    use QuantityTemplates, only: QuantityTemplate_T
    use ReadAPriori, only: read_apriori
    use RetrievalModule, only: Retrieve
    use SpectroscopyCatalog_m, only: Destroy_Line_Database, &
      & Destroy_SpectCat_Database, Spectroscopy
    ! use Test_Parse_Signals_m, only: Test_Parse_Signals
    use Time_M, only: Time_Now
    use Toggles, only: GEN, LEVELS, SWITCHES, TOGGLE
    use Trace_m, only: DEPTH, TRACE_BEGIN, TRACE_END
    use Tree, only: DECORATION, NSONS, SUBTREE
    use VectorsModule, only: DestroyVectorDatabase, DUMP_VECTORS, &
      & Vector_T, VectorTemplate_T
    use VGridsDatabase, only: DestroyVGridDatabase, VGrid_T

    integer, intent(in) ::     ROOT         ! Root of the abstract syntax tree
    integer, intent(out) ::    ERROR_FLAG  ! Nonzero means failure
    integer, intent(in) ::     FIRST_SECTION! Index of son of root of first n_cf
    logical, intent(in) ::     COUNTCHUNKS ! Just count the chunks, print them out and quit
    integer, intent(in) ::     SINGLECHUNK ! Just run this one chunk (0 if all)
    integer, intent(in) ::     LASTCHUNKIN ! Just run range [single,last]
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE

    integer ::                                  chunkNo ! Index of Chunks
    type (MLSChunk_T), dimension(:), pointer :: Chunks  ! of data
    type (FGrid_T), dimension(:), pointer ::     FGrids
    ! Forward model configurations:
    integer ::                                   FIRSTCHUNK ! For chunk loop
    type (ForwardModelConfig_T), dimension(:), &
      & pointer ::                               ForwardModelConfigDatabase
    type (GriddedData_T), dimension(:), &
      & pointer ::                               GriddedDataBase
    type (HGrid_T), dimension(:), pointer ::     HGrids
    integer ::                                   HOWMANY  ! Nsons(Root)
    integer ::                                   I, J     ! Loop inductors
    integer ::                                   LASTCHUNK ! For chunk loop
    type (L1BInfo_T) ::                          L1BInfo  ! File handles etc. for L1B dataset
    type (L2AUXData_T), dimension(:), pointer :: L2AUXDatabase
    type (DirectData_T), dimension(:), pointer :: DirectDatabase
    type (L2GPData_T), dimension(:), pointer  :: L2GPDatabase
    ! type (PCFData_T) ::                          L2pcf
    type (Matrix_Database_T), dimension(:), &
      & pointer ::                               Matrices
    type (TAI93_Range_T) ::                      ProcessingRange  ! Data processing range
    integer ::                                   SON              ! Son of Root
    real    ::                                   t1, t2
    integer ::                                   totalNGC   ! Total num garbage colls.
    integer ::                                   fwmIndex  ! Index
    logical ::                                   show_totalNGC = .true.
    type (Vector_T), dimension(:), pointer ::    Vectors
    type (VGrid_T), dimension(:), pointer ::     VGrids

    ! Arguments for Construct not declared above:
    type (QuantityTemplate_T), dimension(:), pointer :: MifGeolocation
    type (QuantityTemplate_T), dimension(:), pointer :: QtyTemplates
    type (VectorTemplate_T), dimension(:), pointer :: VectorTemplates

    nullify ( chunks, forwardModelConfigDatabase, griddedDataBase, &
      & directDatabase, hGrids, l2auxDatabase, l2gpDatabase, matrices, mifGeolocation, &
      & qtyTemplates, vectors, vectorTemplates, fGrids, vGrids )

    depth = 0
    totalNGC = 0
    if ( toggle(gen) ) call trace_begin ( 'WALK_TREE_TO_DO_MLS_L2', &
      & subtree(first_section,root) )
    call time_now ( t1 )
    call OpenAndInitialize ( processingRange, l1bInfo )
    call add_to_section_timing ( 'open_init', t1)
    i = first_section
    howmany = nsons(root)

    ! -------------------------------------------------------------
    ! ------------------------------------------------------------- Loop over tree

    ! Now loop over the sections in the tree
    do while ( i <= howmany )
      son = subtree(i,root)

      ! First those involved in 'preprocessing'
      select case ( decoration(subtree(1,son)) ) ! section index

        ! --------------------------------------------------------- Init sections
      case ( z_globalsettings )
        call set_global_settings ( son, forwardModelConfigDatabase, fGrids, vGrids, &
          & l2gpDatabase, DirectDatabase, processingRange, l1bInfo )
        call add_to_section_timing ( 'global_settings', t1)
      case ( z_mlsSignals )
        call MLSSignals ( son )
        if ( index(switches,'tps') /= 0 ) then
            ! call test_parse_signals
            call MLSMessage ( MLSMSG_Info, ModuleName, &
              & 'Go back and uncomment the previous line in tree_walker' )
        endif
        call add_to_section_timing ( 'signals', t1)
      case ( z_spectroscopy )
        call spectroscopy ( son )
        call add_to_section_timing ( 'spectroscopy', t1)
      case ( z_readapriori )
        call read_apriori ( son , l2gpDatabase, l2auxDatabase, griddedDataBase )
        call add_to_section_timing ( 'read_apriori', t1)
      case ( z_mergeGrids )
        call mergeGrids ( son, griddedDataBase )

        ! --------------------------------------------------------- Chunk divide
        ! Chunk divide can be a special one, in slave mode, we just listen out
        ! for instructions.
      case ( z_chunkdivide )
        if ( parallel%slave .and. .not. parallel%fwmParallel ) then
          call GetChunkInfoFromMaster ( chunks, chunkNo )
          firstChunk = chunkNo
          lastChunk = chunkNo
          parallel%ChunkNo = chunkNo
        else
          if ( .not. checkPaths) then
            call ChunkDivide ( son, processingRange, l1bInfo, chunks )
            call ComputeAllHGridOffsets ( root, i+1, chunks, l1bInfo, &
            & l2gpDatabase, processingRange )
          else
            allocate(chunks(1), stat=error_flag)
            if ( error_flag /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'unable to allocate chunks' )
          endif
          if ( singleChunk /= 0 ) then
            if ( singleChunk < 0 ) then
              call output ( " single chunk " )
              call output ( singleChunk, advance='yes' )
              call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'single chunk number < 0' )
            else if ( singleChunk > size(chunks) ) then
              call output ( " single chunk " )
              call output ( singleChunk, advance='yes' )
              call output ( " lastChunk " )
              call output ( lastChunk, advance='yes' )
              call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'single chunk number > lastChunk' )
            endif
            firstChunk = singleChunk
            lastChunk = singleChunk
            if ( lastChunkIn > firstChunk ) &
              & lastChunk = min(lastChunkIn, size(chunks))
          else
            firstChunk = 1
            lastChunk = size(chunks)
          end if
          if ( countChunks ) then
            error_flag = 0
            call output ( size(chunks) )
            return
          end if
        end if
        if ( toggle(gen) .and. levels(gen) > 0 ) call dump ( chunks )
        call add_to_section_timing ( 'chunk_divide', t1)

        ! --------------------------------------------------------- Chunk processing
        ! Now construct, fill, join and retrieve live inside the 'chunk loop'
      case ( z_construct, z_fill, z_join, z_retrieve )
        ! Do special stuff in some parallel cases, or where there are
        ! no chunks.
        if ( ( size(chunks) < 1 ) .or. &
          & ( parallel%master .and. .not. parallel%fwmParallel ) .or. &
          & ( parallel%slave .and. parallel%fwmParallel ) ) then
          if ( parallel%master .and. .not. parallel%fwmParallel ) then
            if ( singleChunk /= 0 ) then
              call ReduceChunkDatabase(chunks, singleChunk, lastChunk )
            endif
            call L2MasterTask ( chunks, l2gpDatabase, l2auxDatabase )
          end if
          if ( parallel%slave .and. parallel%fwmParallel ) then
            call ConstructMIFGeolocation ( mifGeolocation, l1bInfo, &
              & chunks(singleChunk) ) 
            call L2FWMSlaveTask ( mifGeolocation )
          end if
          ! Sort out the timings
          select case ( decoration(subtree(1,son)) ) ! section index
          case ( z_construct )
            call add_to_section_timing ( 'construct', t1)
          case ( z_fill )
            call add_to_section_timing ( 'fill', t1)
          case ( z_join )
            call add_to_section_timing ( 'join', t1)
          case ( z_retrieve )
            call add_to_section_timing ( 'retrieve', t1)
          end select
          ! print the timing for FullForwardModel, the following return
	  ! if fmt1 or fmt2 is true
          if ( index(switches, 'fmt') /= 0 .and. &
	     & associated(forwardModelConfigDatabase)) then
                  call printForwardModelTiming ( forwardModelConfigDatabase )
          end if

        else
        ! Otherwise, this is the 'standard' work for these sections.
          do chunkNo = firstChunk, lastChunk ! ----------------------- Chunk loop
            if ( index(switches,'chu') /= 0 ) then
              call output ( " ================ Starting processing for chunk " )
              call output ( chunkNo )
              call output ( " ================ ", advance='yes' )
            end if
            if ( parallel%master .and. parallel%fwmParallel ) &
              & call LaunchFWMSlaves ( chunks ( chunkNo ) )
            j = i
subtrees:   do while ( j <= howmany )
              son = subtree(j,root)
              select case ( decoration(subtree(1,son)) ) ! section index
              case ( z_construct )
                if ( .not. checkPaths) &
                & call MLSL2Construct ( son, l1bInfo, processingRange, &
                  & chunks(chunkNo), qtyTemplates, vectorTemplates, &
                  & fGrids, vGrids, hGrids, l2gpDatabase, forwardModelConfigDatabase, &
                  & mifGeolocation )
                call add_to_section_timing ( 'construct', t1)
              case ( z_fill )
                if ( .not. checkPaths) &
                & call MLSL2Fill ( son, l1bInfo, griddedDataBase, &
                  & vectorTemplates, vectors, qtyTemplates, matrices, vGrids, &
                  & l2gpDatabase, l2auxDatabase, chunks, chunkNo )
                call add_to_section_timing ( 'fill', t1)
              case ( z_join )
                call MLSL2Join ( son, vectors, l2gpDatabase, &
                  & l2auxDatabase, DirectDatabase, chunkNo, chunks )
                call add_to_section_timing ( 'join', t1)
              case ( z_retrieve )
                if ( .not. checkPaths) &
                & call retrieve ( son, vectors, matrices, forwardModelConfigDatabase, &
                  & chunks(chunkNo) )
                call add_to_section_timing ( 'retrieve', t1)
              case default
                exit subtrees
              end select
              j = j + 1
            end do subtrees

            if ( index(switches,'chi1') /= 0 .and. chunkNo > 1) then
              ! Dumps nothing after 1st chunk
            else if ( index(switches,'chi') /= 0 ) then
              call output('Here are our diagnostics for chunk ', advance='no')
              call output(chunkNo, advance='yes')
              call dump_vectors( vectors, details=1, &
              & quantityTypes = (/l_chisqchan, l_chisqmmaf, l_chisqmmif/) )
            end if

            ! Now, if we're dealing with more than one chunk destroy stuff
            ! Otherwise, we'll save them as we may need to output them as l2pc files.
            if ( size(chunks) > 1 .and. &
              & ( singleChunk == 0 .or. lastChunkIn /= 0 ) ) then
              call MLSL2DeConstruct ( qtyTemplates, vectorTemplates, &
                & mifGeolocation, hGrids )
              call DestroyVectorDatabase ( vectors )
              call DestroyMatrixDatabase ( matrices )
              ! if (garbage_collection_by_chunk) call mls_gc_now
              ! if ( index(switches,'ngc') /= 0 ) &
              !  & totalNGC = Say_num_gcs()
            end if
            call ForgetOptimumLon0
            ! print the timing for FullForwardModel
            ! fmt2: at each chunk, fmt1: at last chunk
            if ( index(switches, 'fmt2') /= 0 ) then
                  call printForwardModelTiming ( forwardModelConfigDatabase )
            end if
            if ( index(switches, 'fmt1') /= 0 .and. chunkNo == lastChunk) then
                  call printForwardModelTiming ( forwardModelConfigDatabase )
            end if
            call StripForwardModelConfigDatabase ( forwardModelConfigDatabase )
          end do ! ---------------------------------- End of chunk loop
          ! Clear any locked l2pc bins out.
          call FlushLockedBins
          i = j - 1 ! one gets added back in at the end of the outer loop
        end if

        ! ------------------------------------------- Output section
      case ( z_output ) ! Write out the data
        if ( .not. parallel%slave ) then
          call Output_Close ( son, l2gpDatabase, l2auxDatabase, DirectDatabase, &
            & matrices, size(chunks)==1 .or. singleChunk /= 0 )
        end if

        ! For case where there was one chunk, destroy vectors etc.
        ! This is to guard against destroying stuff needed by l2pc writing
        if ( size(chunks) == 1 .or. &
          & (singleChunk /= 0 .and. lastChunk == 0) ) then
          call MLSL2DeConstruct ( qtyTemplates, vectorTemplates, &
            & mifGeolocation, hGrids )
          call DestroyVectorDatabase ( vectors )
          call DestroyMatrixDatabase ( matrices )
        end if

        ! Now tidy up any remaining `pointer' data.
        ! processingRange needs no deallocation
        if ( index(switches,'gridd') /= 0 .and. .not. parallel%slave &
         & .and. associated(griddedDataBase) ) then
          call Dump(griddedDataBase)
        end if
        call DestroyGriddedDataDatabase ( griddedDataBase )
        call DestroyChunkDatabase (chunks )
        if ( index(switches,'l2gp') /= 0 .and. .not. parallel%slave) then
          call Dump(l2gpDatabase)
        elseif ( index(switches,'cab') /= 0 .and. .not. parallel%slave) then
          call Dump(l2gpDatabase, ColumnsOnly=.true.)
        end if
        call DestroyL2GPDatabase ( l2gpDatabase )
        if ( index(switches,'l2aux') /= 0 .and. .not. parallel%slave) then
          call Dump(l2auxDatabase)
        end if
        call DestroyL2AUXDatabase ( l2auxDatabase )
        call DestroyDirectDatabase ( DirectDatabase )
        ! vectors, vectorTemplates and qtyTemplates destroyed at the
        ! end of each chunk
        call add_to_section_timing ( 'output', t1)

      end select
      i = i + 1
    end do

    ! Now finish up
    call CloseParallel(size(Chunks))
    call destroy_ant_patterns_database
    call DestroyBinSelectorDatabase
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
    call destroyFGridDatabase ( fGrids )
    call DestroyL1BInfo ( l1bInfo )
    error_flag = 0
    if ( toggle(gen) ) call trace_end ( 'WALK_TREE_TO_DO_MLS_L2' )

  contains
    subroutine SayTime ( What )
      character(len=*), intent(in) :: What
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      end if
      call output ( "Timing for " // what // " = " )
      call output ( dble(t2 - t1), advance = 'yes' )
    end subroutine SayTime

    integer function Say_num_gcs ( )
      Say_num_gcs = mls_howmany_gc()
      if ( show_totalNGC ) then
        call output ( "Total = " )
        call output ( Say_num_gcs, advance = 'no' )
        call blanks ( 4, advance = 'no' )
      end if
      call output ( "garbage collections for chunk ")
      call blanks ( 2, advance = 'no' )
      call output ( chunkNo )
      call output ( " = ")
      call output ( Say_num_gcs - totalNGC, advance = 'yes' )
    end function Say_num_gcs

  end subroutine WALK_TREE_TO_DO_MLS_L2

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module TREE_WALKER

! $Log$
! Revision 2.120  2003/12/05 00:42:05  pwagner
! Removed last vestige of explicit garbage collection
!
! Revision 2.119  2003/11/07 00:46:51  pwagner
! New quicker preflight option: --checkPaths
!
! Revision 2.118  2003/11/05 21:27:54  pwagner
! Can enter range of chunks to be processed instead of single
!
! Revision 2.117  2003/10/09 23:31:26  pwagner
! Changed grid switch to gridd
!
! Revision 2.116  2003/09/03 15:56:18  cvuu
! Move the do loop over the forwardModel inside subroutine PrintForwardModelTiming
!
! Revision 2.115  2003/09/02 18:04:38  pwagner
! Can do a singleChunk even if master task
!
! Revision 2.114  2003/08/21 21:24:39  cvuu
! Change the output format for fullForwardModel Timing
!
! Revision 2.113  2003/08/21 16:07:34  livesey
! Now calls FlushLockedBins automatically at the end of each chunk.
!
! Revision 2.112  2003/08/11 23:24:02  pwagner
! Stores ChunkNo as component of L2ParallelInfo_T
!
! Revision 2.111  2003/07/07 23:51:06  pwagner
! Need not pass around l2pc as L2pcf now a saved variable in WriteMetaData
!
! Revision 2.110  2003/06/30 22:56:16  cvuu
! Print mean, std dev for fullForwardModel timing
!
! Revision 2.109  2003/06/24 23:54:07  pwagner
! New db indexes stored for entire direct file
!
! Revision 2.108  2003/06/24 23:00:53  livesey
! Nullified directDatabase
!
! Revision 2.107  2003/06/23 23:55:17  pwagner
! Added DirectData_T to keep track of data written directly
!
! Revision 2.106  2003/06/20 19:38:26  pwagner
! Allows direct writing of output products
!
! Revision 2.105  2003/06/09 22:51:35  pwagner
! Renamed scan_divide to chunk_divide in timings table
!
! Revision 2.104  2003/03/08 00:46:08  pwagner
! Per njl, checks for illegal single chunk runs
!
! Revision 2.103  2003/02/21 21:03:39  pwagner
! Disabled tps switch to eliminate need for Test_Parse_Signals
!
! Revision 2.102  2003/01/13 20:59:14  livesey
! Removed a print statement
!
! Revision 2.101  2002/12/05 19:45:20  pwagner
! Moved MLSFile_T from MLSFiles to MLSCommon
!
! Revision 2.100  2002/12/04 01:18:21  pwagner
! First halting steps toward using filedatabase
!
! Revision 2.99  2002/11/22 12:24:25  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.98  2002/11/22 01:14:06  vsnyder
! Remove USE'd but unreferenced symbols and two unused local variables
!
! Revision 2.97  2002/11/21 18:45:04  livesey
! Changes to way the destroy stuff is called based on chunks
!
! Revision 2.96  2002/11/21 18:38:34  livesey
! Bug fix in canWriteL2PC flag passed to output_close
!
! Revision 2.95  2002/10/08 17:41:50  livesey
! Various bug fixes associated with FWMParallel
!
! Revision 2.94  2002/10/08 17:36:23  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.93  2002/10/05 00:44:29  livesey
! Included the FWMParallel stuff
!
! Revision 2.92  2002/09/25 20:09:32  livesey
! Changes to allow chunk based forwardModel configs
!
! Revision 2.91  2002/09/24 22:17:50  pwagner
! Defines t1 before call to add_to_section_timing
!
! Revision 2.90  2002/09/24 18:18:47  pwagner
! Consistent with add_to_section_timing now calling time_now at its end
!
! Revision 2.89  2002/08/28 22:29:19  pwagner
! Moved DestroyL1BInfo to after end of loop of chunks
!
! Revision 2.88  2002/08/22 01:26:00  vsnyder
! Cosmetic changes
!
! Revision 2.87  2002/08/21 19:05:04  livesey
! Removed calls to H5Open/close as they are now done in MLSL2.
!
! Revision 2.86  2002/08/21 00:55:35  livesey
! Changed error to warning when fail to close hdf5 (at least for the
! moment).
!
! Revision 2.85  2002/08/20 22:10:50  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.84  2002/08/07 00:05:27  livesey
! Added calls to H5Open_F and H5Close_F
!
! Revision 2.83  2002/08/04 16:10:43  mjf
! Added some nullify statements for Sun's rubbish compiler.
!
! Revision 2.82  2002/03/20 00:49:44  pwagner
! chi1 dumps only 1st chunk, unlike chi
!
! Revision 2.81  2002/02/05 00:44:03  pwagner
! Added garbage collection stuff
!
! Revision 2.80  2002/01/24 00:58:28  livesey
! Now calls MergeGrids at the appropriate time
!
! Revision 2.79  2002/01/22 18:14:47  livesey
! Fixed typo
!
! Revision 2.78  2002/01/21 23:11:06  livesey
! Added call to DestroyBinSelectorsDatabase etc.
!
! Revision 2.77  2002/01/18 18:55:37  livesey
! Code to support the --chunk option
!
! Revision 2.76  2002/01/09 22:56:17  livesey
! Now sends slaves all chunks as regular HGrids need them.
!
! Revision 2.75  2001/12/16 00:57:12  livesey
! Now passes all chunks to construct
!
! Revision 2.74  2001/12/14 01:43:20  livesey
! Passes processingRange to HGrid via Construct
!
! Revision 2.73  2001/12/13 23:21:26  livesey
! Added countChunks option
!
! Revision 2.72  2001/12/10 20:22:09  livesey
! Added code for EmpiricalGeometry.
!
! Revision 2.71  2001/11/20 00:48:15  livesey
! Fixed problem when no chunks to process
!
! Revision 2.70  2001/11/16 17:24:13  livesey
! Now calls new ChunkDivide routine.
!
! Revision 2.69  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.68  2001/10/31 19:07:52  livesey
! Hooked fGrids into quantity templates
!
! Revision 2.67  2001/10/31 18:36:58  livesey
! Added fGrids stuff
!
! Revision 2.66  2001/10/26 23:16:15  pwagner
! Similar dump interfaces for l2gp, l2aux, Griddeddata databases
!
! Revision 2.65  2001/10/12 23:14:22  pwagner
! Debugging when dumping diagnostics; may remove later
!
! Revision 2.64  2001/10/09 23:43:42  pwagner
! Some further improvements in dumping vectors
!
! Revision 2.63  2001/10/02 16:49:56  livesey
! Removed fmStat%finished and change loop ordering in forward models
!
! Revision 2.62  2001/09/29 00:01:00  pwagner
! Fixed various timing problems
!
! Revision 2.61  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
! Revision 2.60  2001/09/21 17:40:57  pwagner
! May dump diagnostic quantities chi..
!
! Revision 2.59  2001/09/10 23:37:32  livesey
! New GriddedData etc.
!
! Revision 2.58  2001/08/06 18:34:59  pwagner
! Now dumps l2gp and l2aux databases when asked
!
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
