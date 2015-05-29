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
MODULE Construct                ! The construct module for the MLS L2 sw.
!=============================================================================

  ! This module performs the `construct' task for the level 2 software.  This
  ! task involves constructing templates for vector quantities, vectors and
  ! matrices.

  implicit none

  private

  public :: MLSL2Construct, MLSL2DeConstruct, ConstructMIFGeolocation
  
  !------------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName="$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! --------------------------------------------- ConstructMIFGeolocation --
  subroutine ConstructMIFGeolocation ( mifGeolocation, filedatabase, chunk )
    ! mifGeolocation is just quantity templates containing geolocation
    ! information for the GHz and THz modules.  The software can then
    ! point to these for geolocation information for all minor frame
    ! quantities saving file IO and memory.
    use CHUNKS_M, only: MLSCHUNK_T
    use CONSTRUCTQUANTITYTEMPLATES, only: CONSTRUCTMINORFRAMEQUANTITY
    use QUANTITYTEMPLATES, only: QUANTITYTEMPLATE_T
    use MLSCOMMON, only: MLSFILE_T
    use MLSSIGNALS_M, only: MODULES
    use MLSMESSAGEMODULE, only: MLSMSG_ERROR, MLSMSG_ALLOCATE, &
      & MLSMESSAGE

    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type (MLSChunk_T), intent(in) :: chunk
    
    ! Local variables
    integer :: INSTRUMENTMODULEINDEX    ! Loop counter
    integer :: STATUS                   ! Flag

    if ( .not. associated ( mifGeolocation ) ) then
      ! Don't overwrite it if we already have it, e.g. from previous construct
      ! or forge.
      allocate ( mifGeolocation(size(modules)), STAT=status )
      if ( status/=0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//"mifGeolocation" )
    
      ! Now try to fill it if we have any L1BFiles
      if (associated(filedatabase) ) then
        do instrumentModuleIndex = 1, size(modules)
          call ConstructMinorFrameQuantity ( instrumentModuleIndex, &
          mifGeolocation(instrumentModuleIndex), &
          filedatabase=filedatabase, chunk=chunk)
        end do
      else
        mifGeolocation%noSurfs = 0
        mifGeolocation%noInstances = 0
      end if
    end if
  end subroutine ConstructMIFGeolocation

  ! ---------------------------------------------  MLSL2Construct  -----
  subroutine MLSL2Construct ( root, filedatabase, processingRange, chunk, &
       & quantityTemplatesBase, vectorTemplates, vectors, FGrids, HGrids, &
       & l2gpDatabase, ForwardModelConfigDatabase, griddedDataBase, &
       & mifGeolocation )

  ! This is the `main' subroutine for this module

    use CHUNKS_M, only: MLSCHUNK_T
    use CONSTRUCTQUANTITYTEMPLATES, only: &
      & CREATEQTYTEMPLATEFROMMLSCFINFO, FORGEMINORFRAMES
    use CONSTRUCTVECTORTEMPLATES, only: CREATEVECTEMPLATEFROMMLSCFINFO
    use DUMPCOMMAND_M, only: BOOLEANFROMANYGOODRADIANCES, &
      & BOOLEANFROMANYGOODVALUES, &
      & BOOLEANFROMCATCHWARNING, BOOLEANFROMCOMPARINGQTYS, BOOLEANFROMFORMULA, &
      & DUMPCOMMAND
    use FGRID, only: FGRID_T
    use FORWARDMODELCONFIG, only: ADDFORWARDMODELCONFIGTODATABASE, &
      & FORWARDMODELCONFIG_T
    use FORWARDMODELSUPPORT, only: CONSTRUCTFORWARDMODELCONFIG
    use GRIDDEDDATA, only: GRIDDEDDATA_T
    use HGRIDSDATABASE, only: ADDHGRIDTODATABASE, HGRID_T
    use HGRID, only: CREATEHGRIDFROMMLSCFINFO
    use INIT_TABLES_MODULE, only: S_ANYGOODVALUES, S_ANYGOODRADIANCES, &
      & S_BOOLEAN, S_CATCHWARNING, S_COMPARE, S_DUMP, &
      & S_FORGE, S_FORWARDMODEL, S_HGRID, &
      & S_PHASE, S_QUANTITY, S_REEVALUATE, S_TIME, S_VECTORTEMPLATE
    use L2GPDATA, only: L2GPDATA_T
    use MLSCOMMON, only: MLSFILE_T, TAI93_RANGE_T
    use MLSL2OPTIONS, only: L2CFNODE, NEED_L1BFILES, SPECIALDUMPFILE
    use MLSL2TIMINGS, only: SECTION_TIMES, TOTAL_TIMES, ADDPHASETOPHASENAMES
    use MLSMESSAGEMODULE, only: MLSMESSAGERESET
    use MORETREE, only: Get_Label_And_Spec, GET_SPEC_ID
    use Next_Tree_Node_m, only: Next_Tree_Node, Next_Tree_Node_State
    use OUTPUT_M, only: BLANKS, OUTPUT, &
      & REVERTOUTPUT, SWITCHOUTPUT
    use QUANTITYTEMPLATES, only: ADDQUANTITYTEMPLATETODATABASE, &
      & QUANTITYTEMPLATE_T
    use TIME_M, only: TIME_NOW
    use TOGGLES, only: GEN, LEVELS, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATE
    use VECTORSMODULE, only: ADDVECTORTEMPLATETODATABASE, &
      & VECTOR_T, VECTORTEMPLATE_T

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Root of the tree for the Construct section
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type (TAI93_Range_T), intent(in) :: processingRange
    type (MLSChunk_T), intent(in) :: chunk
    type (QuantityTemplate_T), dimension(:), pointer :: quantityTemplatesBase
    type (VectorTemplate_T), dimension(:), pointer :: vectorTemplates
    type (vector_T), dimension(:), pointer :: Vectors
    type (FGrid_T), dimension(:), pointer :: fGrids
    type (HGrid_T), dimension(:), pointer :: HGrids
    type (L2GPData_T), dimension(:), pointer :: L2GPDatabase
    type (ForwardModelConfig_T), dimension(:), pointer :: ForwardModelConfigDatabase
    type (GriddedData_T), dimension(:), pointer        :: griddedDataBase
    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation

    ! Local variables

    integer :: KEY              ! S_... from Init_Tables_Module.
    integer :: Me = -1          ! String index for trace
    integer :: Me_Spec = -1     ! String index for trace
    integer :: NAME             ! Sub-rosa index of name
    integer :: SON              ! Son or grandson of Root
    type(next_tree_node_state) :: State ! of tree traverser
    real :: T1, T2              ! for timing
    logical :: TIMING

    ! Executable code
    timing = section_times
    if ( timing ) call time_now ( t1 )

    call trace_begin ( me, "MLSL2Construct", root, cond=toggle(gen) )
    if ( specialDumpFile /= ' ' ) &
      & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )

    ! First we're going to setup our mifGeolocation quantityTemplates.
    if ( NEED_L1BFILES ) &
      & call ConstructMIFGeolocation ( mifGeolocation, filedatabase, chunk )

    ! The rest is fairly simple really.  We just loop over the mlscf 
    ! instructions and hand them off to people

    do
      son = next_tree_node ( root, state )
      if ( son == 0 ) exit
      call trace_begin ( me_spec, "Construct.spec", son, &
        & cond=toggle(gen) .and. levels(gen) > 0 )
      call get_label_and_spec ( son, name, key )
      L2CFNODE = key

      ! Node_id(key) is now n_spec_args.

      if ( get_spec_id(key) /= s_catchWarning ) &
        & call MLSMessageReset( clearLastWarning=.true. )
      select case( get_spec_id(key) )
      case ( s_anygoodvalues )
        call decorate ( key, &
          & BooleanFromAnyGoodValues ( key, vectors ) )
      case ( s_anygoodradiances )
        call decorate ( key, &
          & BooleanFromAnyGoodRadiances ( key, chunk, filedatabase ) )
      case ( s_catchWarning )
        call decorate ( key,  BooleanFromCatchWarning ( key ) )
      case ( s_compare )
        call decorate ( key,  BooleanFromComparingQtys ( key, vectors ) )
      case ( s_dump )
        call dumpCommand ( key, quantityTemplatesBase, &
          & vectorTemplates, forwardModelConfigs=forwardModelConfigDatabase, &
          & hGrids=hGrids, griddedDataBase=griddedDataBase )
      case ( s_forge )
        call ForgeMinorFrames ( key, mifGeolocation )
      case ( s_forwardModel )
        call decorate ( key, AddForwardModelConfigToDatabase ( &
          & forwardModelConfigDatabase, &
          & ConstructForwardModelConfig ( name, key, .false. ) ) )
      case ( s_hgrid )
        call decorate ( key, AddHGridToDatabase ( hGrids, &
          & CreateHGridFromMLSCFInfo ( name, key, filedatabase, l2gpDatabase, &
          & processingRange, chunk ) ) )
      case ( s_phase )
        call addPhaseToPhaseNames ( name, key )

      case ( s_quantity )
        call decorate ( key, AddQuantityTemplateToDatabase ( &
          & quantityTemplatesBase, CreateQtyTemplateFromMLSCfInfo ( name, key, &
            & fGrids, hGrids, filedatabase, chunk, mifGeolocation ) ) )
      case ( s_Boolean )
        call decorate ( key,  BooleanFromFormula ( name, key ) )
      case ( s_Reevaluate )
        call decorate ( key,  BooleanFromFormula ( 0, key ) )
      case ( s_vectortemplate )
        call decorate ( key, AddVectorTemplateToDatabase ( vectorTemplates, &
          & CreateVecTemplateFromMLSCfInfo ( name, key, quantityTemplatesBase ) ) )
      case ( s_time )
        if ( timing ) then
          call sayTime
        else
          call time_now ( t1 )
          timing = .true.
        end if
      case default ! Can't get here if tree_checker worked correctly
      end select
      call trace_end ( "Construct.spec", cond=toggle(gen) .and. levels(gen) > 0 )
    end do

    if ( specialDumpFile /= ' ' ) call revertOutput
    call trace_end ( "MLSL2Construct", cond=toggle(gen) )

    if ( timing ) call sayTIme

  contains
    subroutine SayTime
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      end if
      call output ( "Timing for MLSL2Construct = " )
      call output ( DBLE(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime
  end subroutine MLSL2Construct

  ! -------------------------------------------  MLSL2DeConstruct  -----
  subroutine MLSL2DeConstruct ( quantityTemplatesBase, vectorTemplates, &
    &                           mifGeolocation, hGrids )

  ! DeConstruct the Quantity and Vector template databases.

    use HGRIDSDATABASE, only: DESTROYHGRIDDATABASE, HGRID_T
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use OUTPUT_M, only: OUTPUT
    use QUANTITYTEMPLATES, only: DESTROYQUANTITYTEMPLATEDATABASE, &
      & QUANTITYTEMPLATE_T
    use TOGGLES, only: SWITCHES
    use VECTORSMODULE, only: DESTROYVECTORTEMPLATEDATABASE, VECTORTEMPLATE_T

    type (QuantityTemplate_T), dimension(:), pointer :: quantityTemplatesBase
    type (VectorTemplate_T), dimension(:), pointer :: vectorTemplates
    type (QuantityTemplate_T), dimension(:), pointer :: mifGeolocation
    type (HGrid_T), dimension(:), pointer :: hGrids
    logical :: verbose

    ! Executable code
    verbose = ( switchDetail(switches, 'destroy' ) > -1 )

    if ( verbose ) call output( 'About to destroy vectortemplate db', advance='yes' )
    call destroyVectorTemplateDatabase ( vectorTemplates )
    if ( verbose ) call output( 'About to destroy quantitytemplate db', advance='yes' )
    call destroyQuantityTemplateDatabase ( quantityTemplatesBase )
    if ( verbose ) call output( 'About to destroy mifGeolocation db', advance='yes' )
    call destroyQuantityTemplateDatabase ( mifGeolocation )
    if ( verbose ) call output( 'About to destroy hGrid db', advance='yes' )
    call destroyHGridDatabase ( hGrids )
  end subroutine MLSL2DeConstruct

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

END MODULE Construct
!=============================================================================

!
! $Log$
! Revision 2.73  2015/05/29 17:48:44  vsnyder
! Remove 'chunk' argument from MIFGeolocation
!
! Revision 2.72  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.71  2013/12/12 02:11:26  vsnyder
! Use iterator to handle variables, and IF and SELECT constructs
!
! Revision 2.70  2013/10/09 23:41:00  vsnyder
! Add Evaluate_Variable
!
! Revision 2.69  2013/08/30 02:45:35  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.68  2013/08/21 00:23:01  pwagner
! -g[level] can trace individual Construct specs
!
! Revision 2.67  2012/08/16 17:55:00  pwagner
! Exploit level 2-savvy MLSMessage
!
! Revision 2.66  2012/01/05 01:20:17  pwagner
! Capitalized USEd stuff
!
! Revision 2.65  2011/06/29 21:54:51  pwagner
! Some cases may safely omit l1b files
!
! Revision 2.64  2010/04/05 17:32:04  honghanh
! Make filedatabase and chunk in ConstructMinorFrameQuantity optional
!
! Revision 2.63  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.62  2007/12/07 01:13:18  pwagner
! Lets us catch warnings and assign to runtime Booleans
!
! Revision 2.61  2007/11/15 22:51:10  pwagner
! Boolean functions moved to DumpCommand
!
! Revision 2.60  2007/10/24 00:14:58  pwagner
! Removed unused declarations
!
! Revision 2.59  2007/03/23 00:24:12  pwagner
! Switch destroy warns when destroying dbs
!
! Revision 2.58  2006/06/12 16:28:56  pwagner
! Added ability to dump Gridded Data
!
! Revision 2.57  2006/03/07 00:51:32  pwagner
! May change already-set Booleans via reevaluate command
!
! Revision 2.56  2006/03/04 00:20:45  pwagner
! May skip retrieval, directWrites depending on runtime Booleans
!
! Revision 2.55  2006/02/10 21:19:07  pwagner
! dumps may go to special dumpfile
!
! Revision 2.54  2006/01/06 01:16:34  pwagner
! silent boolean field can silence selected phases
!
! Revision 2.53  2006/01/04 01:27:11  vsnyder
! Comment out use for unreference L1BInfo_T
!
! Revision 2.52  2005/06/03 02:05:29  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! get VGrids from VGridsDatabase instead of passing as an argument.
!
! Revision 2.51  2005/05/31 17:51:16  pwagner
! Began switch from passing file handles to passing MLSFiles
!
! Revision 2.50  2005/03/12 00:50:27  pwagner
! May restart warnings counter at each phase
!
! Revision 2.49  2004/10/13 02:24:02  livesey
! Added vGrids to Forge command
!
! Revision 2.48  2004/05/20 19:47:55  vsnyder
! Do all dumping by way of DumpCommand
!
! Revision 2.47  2004/05/19 19:16:09  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.46  2004/05/18 01:24:31  vsnyder
! Add HGrids argument to DumpCommand
!
! Revision 2.45  2004/05/11 02:54:36  vsnyder
! Remove USEs and declarations for unreferenced symbols
!
! Revision 2.44  2004/05/03 16:37:11  pwagner
! Get dump from QuantityTemplates module
!
! Revision 2.43  2004/05/01 04:05:26  vsnyder
! Add Dump command
!
! Revision 2.42  2003/10/22 21:17:06  pwagner
! aPhaseName: Phase added to Fill, Construct sections to time phases
!
! Revision 2.41  2003/07/15 18:18:16  livesey
! Change to forward model config call
!
! Revision 2.40  2003/06/20 19:37:06  pwagner
! Quanities now share grids stored separately in databses
!
! Revision 2.39  2003/05/28 04:39:32  livesey
! Removed some obsolete checking
!
! Revision 2.38  2002/10/08 17:36:19  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.37  2002/10/05 00:43:47  livesey
! Split out ConstructMIFGeolocation
!
! Revision 2.36  2002/09/26 23:58:57  livesey
! Bug fix, decorated wrong node
!
! Revision 2.35  2002/09/25 20:07:41  livesey
! Can now construct forward models inside construct
!
! Revision 2.34  2002/09/18 22:49:42  pwagner
! Receives returnStatus from CreateQtyTemplateFromMLSCFInfo
!
! Revision 2.33  2002/08/20 22:43:37  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.32  2001/12/16 00:57:00  livesey
! Now takes all chunks as argument as HGrid needs them
!
! Revision 2.31  2001/12/14 01:42:47  livesey
! Passes processingRange to HGrid construction
!
! Revision 2.30  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.29  2001/10/31 19:07:15  livesey
! Added fGrid stuff
!
! Revision 2.28  2001/10/05 23:30:37  pwagner
! Fixed small bug in destroyingQuantTempldb
!
! Revision 2.27  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.26  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
! Revision 2.25  2001/05/21 20:57:55  livesey
! Fixed bug, was overwriting mifGelocation with each new construct.
!
! Revision 2.24  2001/05/12 00:16:55  livesey
! Big fix, only dump hGrids etc. if they exist.
!
! Revision 2.23  2001/05/10 00:43:12  livesey
! Moved ownership of hGrids into tree walker to allow multiple calls
!
! Revision 2.22  2001/04/28 01:43:10  vsnyder
! Improved the timing message
!
! Revision 2.21  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.20  2001/04/25 19:29:31  livesey
! Fixed bug with forge, now sets mafIndex and mafCounter appropriately
!
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

