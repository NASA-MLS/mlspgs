! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Join                     ! Join together chunk based data.
!=============================================================================

  ! This module performs the 'join' task in the MLS level 2 software.

  implicit none
  private
  public :: MLSL2Join, JoinL2GPQuantities, JoinL2AuxQuantities

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! logical, parameter, private :: DEEBUG = .true.           ! Usually FALSE
  logical, parameter, private :: FORCEDIRWRITEREOPEN = .false. ! Usually FALSE
  logical, parameter, private :: SKIPMETADATA = .false. ! Usually FALSE

  ! Parameters for Announce_Error

  integer :: ERROR
  integer, parameter :: NO_ERROR_CODE=0
  integer, parameter :: NotAllowed=1

contains ! =====     Public Procedures     =============================

  ! --------------------------------------------------  MLSL2Join  -----

  ! This is the main routine for the Join block.  This one just goes
  ! through the tree and dispatches work to other routines.

  subroutine MLSL2Join ( root, vectors, l2gpDatabase, l2auxDatabase, &
    & DirectDataBase, chunkNo, chunks )
    ! Imports
    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use DirectWrite_m, only: DirectData_T
    use Init_Tables_Module, only: S_L2GP, S_L2AUX, S_TIME, S_DIRECTWRITE, S_LABEL
    use L2GPData, only: L2GPDATA_T
    use L2AUXData, only: L2AUXDATA_T
    use L2ParInfo, only: PARALLEL, WAITFORDIRECTWRITEPERMISSION
    use MLSCommon, only: MLSCHUNK_T
    use MLSL2Options, only: CHECKPATHS
    use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES, &
      & add_to_directwrite_timing, add_to_section_timing
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MoreTree, only: GET_SPEC_ID
    use Output_m, only: OUTPUT, BLANKS
    use TOGGLES, only: GEN, TOGGLE, LEVELS, SWITCHES
    use Tree, only: SUBTREE, NSONS, NODE_ID
    use TREE_TYPES, only: N_NAMED
    use Time_M, only: Time_Now
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use VectorsModule, only: VECTOR_T
    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the JOIN section in the AST
    type (Vector_T), dimension(:), pointer :: vectors
    type (L2GPData_T), dimension(:), pointer :: l2gpDatabase
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase
    type (DirectData_T), dimension(:), pointer :: DirectDatabase
    integer, intent(in) :: chunkNo
    type (MLSChunk_T), dimension(:), intent(in) :: chunks

    ! Local parameters
    integer, parameter :: DELAY = 500000  ! For Usleep, no. microsecs
    ! External (C) function
    external :: Usleep
    ! Local variables
    integer :: DIRECTWRITENODEGRANTED   ! Which request was granted
    integer :: DWINDEX                  ! Direct Write index
    real :: DWT1                        ! Time we started
    real :: DWT2                        ! Time we finished
    integer :: KEY                      ! Tree node
    integer :: MLSCFLINE                ! Line number in l2cf
    integer :: NODIRECTWRITES           ! Array size
    integer :: NODIRECTWRITESCOMPLETED  ! Counter
    integer :: PASS                     ! Loop counter
    integer :: SON                      ! Tree node
    integer :: SPECID                   ! Type of l2cf line this is
    integer :: TICKET                   ! Direct write permission ticket
    logical :: CREATEFILE               ! Flag
    logical :: TIMING                   ! Flag
    real :: T1                          ! Time we started
    real :: T2                          ! Time we finished

    integer :: I
    logical :: DEEBUG
    
    ! Executable code
    DEEBUG = (index(switches, 'direct') /= 0)
    if ( toggle(gen) ) call trace_begin ( "MLSL2Join", root )
    timing = section_times
    if ( timing ) call time_now ( t1 )

    ! This is going to be somewhat atypical, as the code may run in 'passes'
    ! In pass 1 we do all the regular join statements and count the direct writes
    ! In pass 2 we log all our direct write requests.
    ! In the later passes (as many as there are direct writes) we do the
    ! direct writes we've been given permission for.

    ! In the non-parallel slave mode, one pass is sufficient.

    error = 0
    pass = 1
    noDirectWrites = 0
    noDirectWritesCompleted = 0
    ticket = 0                          ! Default value for serial case
    passLoop: do
      ! For the later passes, we wait for permission to do one of the direct writes
      if ( pass > 2 ) then
        call time_now ( dwt2 )
        call WaitForDirectWritePermission ( directWriteNodeGranted, ticket, createFile )
        call output ( "Got permission for ticket " )
        call output ( ticket )
        call output ( " node " )
        call output ( directWriteNodeGranted, advance='yes' )
        call add_to_directwrite_timing ( 'waiting', dwt2)
      end if
      
      ! Simply loop over lines in the l2cf
      do mlscfLine = 2, nsons(root) - 1 ! Skip begin/end section
        son = subtree(mlscfLine,root)
        if ( node_id(son) == n_named ) then ! Is spec labeled?
          key = subtree(2,son)
        else
          key = son
        end if
        specId = get_spec_id ( key )
        select case ( specId )
        case ( s_time )
          ! Only say the time the first time round
          if ( pass == 1 ) then
            if ( timing ) then
              call sayTime
            else
              call time_now ( t1 )
              timing = .true.
            end if
          end if
        case ( s_l2gp, s_l2aux )
          ! Only do these the first time round
          if ( pass == 1 .and. .not. checkpaths ) then
            call JoinQuantities ( son, vectors, l2gpDatabase, l2auxDatabase, &
              & chunkNo, chunks )
          end if
        case ( s_label )
          ! Only do these the first time round
          if ( pass == 1 .and. .not. checkpaths ) then
            call LabelVectorQuantity ( son, vectors )
          end if
        case ( s_directWrite )
          call time_now ( dwt1 )
          if ( pass == 1 ) then
            ! On the first pass just count the number of direct writes
            noDirectWrites = noDirectWrites + 1
            ! Unless we're not a slave in which case just get on with it.
            if ( .not. parallel%slave ) then
              if(DEEBUG)print*,'Calling DirectWrite, not slave'
              call time_now ( dwt2 )
              call DirectWriteCommand ( son, ticket, vectors, DirectdataBase, &
                & chunkNo, chunks )
              call add_to_directwrite_timing ( 'writing', dwt2)
            end if
          else if ( pass == 2 ) then
            ! On the second pass, log all our direct write requests.
            if(DEEBUG)print*,'Calling direct write to do a setup'
            call DirectWriteCommand ( son, ticket, vectors, DirectdataBase, &
              & chunkNo, chunks, makeRequest=.true. )
          else
            ! On the later passes we do the actual direct write we've been
            ! given permission for.
            if ( son == directWriteNodeGranted ) then
              call time_now ( dwt2 )
              if(DEEBUG)print*,'Calling direct write to do the write'
              call DirectWriteCommand ( son, ticket, vectors, DirectdataBase, &
                & chunkNo, chunks, create=createFile )
              call add_to_directwrite_timing ( 'writing', dwt2)
              noDirectWritesCompleted = noDirectWritesCompleted + 1
              ! If that was the last one then bail out
              if ( noDirectWritesCompleted == noDirectWrites ) exit passLoop
            end if
          end if
          call add_to_section_timing ( 'directwrite', dwt1)
        end select
      end do                            ! End loop over l2cf lines

      ! If we're not in parallel mode then one pass is enough
      if ( .not. parallel%slave ) exit passLoop

      ! Bail out of pass loop if there are no direct writes, or there was
      ! an error.
      if ( noDirectWrites == 0 .or. error /= 0 ) exit passLoop
      pass = pass + 1
    end do passLoop                     ! End loop over passes

    ! Check for errors
    if ( error /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Error in Join section' )

    if ( toggle(gen) ) call trace_end ( "MLSL2Join" )
    if ( timing ) call sayTime

  contains
    ! Private procedure
    subroutine SayTime
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      end if
      call output ( "Timing for MLSL2Join =" )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

  end subroutine MLSL2Join

  ! ------------------------------------------------ DirectWriteCommand -----
  subroutine DirectWriteCommand ( node, ticket, vectors, DirectDataBase, &
    & chunkNo, chunks, makeRequest, create )
    ! Imports
    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use DirectWrite_m, only: DirectData_T, &
      & AddDirectToDatabase, DirectWrite_l2GP, DirectWrite_l2aux, &
      & SetUpNewDirect
    use Expr_m, only: EXPR
    use Hdf, only: DFACC_CREATE, DFACC_RDONLY, DFACC_RDWR
    use Init_tables_module, only: F_SOURCE, F_PRECISION, F_HDFVERSION, F_FILE, F_TYPE
    use Init_tables_module, only: L_L2GP, L_L2AUX, L_L2DGG, L_L2FWM, &
      & L_PRESSURE, L_ZETA
    use intrinsic, only: L_NONE, L_GEODANGLE, L_HDF, L_SWATH, &
      & L_MAF, PHYQ_DIMENSIONLESS
    use L2ParInfo, only: PARALLEL, LOGDIRECTWRITEREQUEST, FINISHEDDIRECTWRITE
    use MLSCommon, only: MLSCHUNK_T, R4, R8, RV
    use MLSFiles, only: HDFVERSION_5, NAMENOTFOUND, &
      & MLS_EXISTS, split_path_name, GetPCFromRef, &
      & mls_io_gen_openF, mls_io_gen_closeF, mls_sfstart, mls_sfend
    use MLSHDFEOS, only: mls_swath_in_file
    use MLSL2Options, only: CHECKPATHS, TOOLKIT, DEFAULT_HDFVERSION_WRITE
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSPCF2, only: mlspcf_l2gp_start, mlspcf_l2gp_end, &
      & mlspcf_l2dgm_start, mlspcf_l2dgm_end, mlspcf_l2fwm_full_start, &
      & mlspcf_l2fwm_full_end, &
      & mlspcf_l2dgg_start, mlspcf_l2dgg_end
    use MoreTree, only: GET_FIELD_ID
    use Output_m, only: OUTPUT, BLANKS
    use OutputAndClose, only: add_metadata
    use String_Table, only: DISPLAY_STRING, GET_STRING
    use TOGGLES, only: GEN, TOGGLE, LEVELS, SWITCHES
    use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, NULL_TREE, SOURCE_REF, &
      & SUB_ROSA, SUBTREE
    use VectorsModule, only: VECTOR_T, VECTORVALUE_T, VALIDATEVECTORQUANTITY, &
      & GETVECTORQTYBYTEMPLATEINDEX
    ! Dummy arguments
    integer, intent(in) :: NODE         ! Of the JOIN section in the AST
    integer, intent(in) :: TICKET       ! Ticket number from master
    type (Vector_T), dimension(:), pointer :: VECTORS
    type (DirectData_T), dimension(:), pointer :: DirectDatabase
    integer, intent(in) :: CHUNKNO
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    logical, intent(in), optional :: MAKEREQUEST  ! Set on first pass through
    logical, intent(in), optional :: CREATE       ! Set if slave is to create file.
    ! Local parameters
    integer, parameter :: MAXFILES = 1000 ! Set for an internal array
    ! Saved variable - used to work out the createFile flag in the serial case
    integer, dimension(maxFiles), save :: CREATEDFILENAMES = 0
    integer, save :: NOCREATEDFILES=0   ! Number of files created
    ! Local variables
    integer :: fileaccess               ! DFACC_CREATE or DFACC_RDWR
    integer :: DBINDEX
    integer :: ERRORTYPE
    integer :: EXPECTEDTYPE             ! l2gp/l2aux
    integer :: FIELDINDEX               ! Type of field in l2cf line
    integer :: FILE                     ! File name string index
    integer :: FILETYPE
    integer :: GSON                     ! Son of son
    integer :: HANDLE                   ! File handle from hdf/hdf-eos
    integer :: HDFNAMEINDEX             ! String index for output name
    integer :: HDFVERSION               ! 4 or 5
    integer :: KEYNO                    ! Loop counter, field in l2cf line
    integer :: LASTFIELDINDEX           ! Type of previous field in l2cf line
    integer :: NOSOURCES                ! No. things to output
    integer :: OUTPUTTYPE               ! l_l2gp, l_l2aux, l_l2fwm, l_l2dgg
    integer :: SON                      ! A tree node
    integer :: SOURCE                   ! Loop counter
    integer :: RETURNSTATUS
    logical :: CREATEFILEFLAG           ! Flag (often copy of create)
    logical :: MYMAKEREQUEST            ! Copy of makeRequest
    logical :: DUMMY
    logical, dimension(:), pointer :: CREATETHISSOURCE
    character(len=256), dimension(:), pointer :: NAMEBUFFER
    integer :: record_length
    integer :: l2gp_Version

    integer :: EXPRUNITS(2)             ! From expr
    real (r8) :: EXPRVALUE(2)           ! From expr
    integer, dimension(:), pointer :: SOURCEVECTORS ! Indicies
    integer, dimension(:), pointer :: SOURCEQUANTITIES ! Indicies
    integer, dimension(:), pointer :: PRECISIONVECTORS ! Indicies
    integer, dimension(:), pointer :: PRECISIONQUANTITIES ! Indicies
    character(len=1024) :: FILENAME     ! Output full filename
    character(len=1024) :: FILE_BASE    ! made up of
    character(len=1024) :: PATH         ! path/file_base
    character(len=1024) :: HDFNAME      ! Output swath/sd name
    type(VectorValue_T), pointer :: QTY ! The quantity
    type(VectorValue_T), pointer :: PRECQTY ! The quantities precision
    type(DirectData_T) :: newDirect
    type(DirectData_T), pointer :: thisDirect
    logical :: DEEBUG
    integer, external :: he5_SWclose

    ! Executable code
    DEEBUG = (index(switches, 'direct') /= 0)

    myMakeRequest = .false.
    if ( present ( makeRequest ) ) myMakeRequest = makeRequest

    ! Direct write is probably going to be come the predominant form
    ! of output in the software, as the other forms have become a
    ! little too intensive.

    ! This routine will be called twice for each direct write.  The first
    ! (makeRequests=.true.) is a pass through to check the syntax etc.
    ! If it is sucessfull a request will be made to the master.

    ! The second call will be to actually do the writing.
    lastFieldIndex = 0
    noSources = 0
    hdfVersion = DEFAULT_HDFVERSION_WRITE
    do keyNo = 2, nsons(node)           ! Skip DirectWrite command
      son = subtree ( keyNo, node )
      if ( keyNo > 2 ) lastFieldIndex = fieldIndex
      fieldIndex = get_field_id ( son )
      select case ( fieldIndex )
      case ( f_source )
        noSources = noSources + 1
      case ( f_precision )
        if ( lastFieldIndex /= f_source ) call Announce_Error ( son, no_error_code, &
            & 'A precision can only be given immediately following a source' )
      case ( f_hdfVersion )
        call expr ( subtree(2,son), exprUnits, exprValue )
        if ( exprUnits(1) /= phyq_dimensionless ) &
          & call Announce_error ( son, NO_ERROR_CODE, &
          & 'No units allowed for hdfVersion: just integer 4 or 5')
        hdfVersion = exprValue(1)
      case ( f_file )
        file = sub_rosa(subtree(2,son))
      case ( f_type )
        outputType = decoration(subtree(2,son))
      end select
    end do

    call get_string ( file, filename, strip=.true. )
    
    ! Now identify the quantities we're after
    nullify ( sourceVectors, sourceQuantities, precisionVectors, precisionQuantities )
    call Allocate_test ( sourceVectors, noSources, 'sourceVectors', ModuleName )
    call Allocate_test ( sourceQuantities, noSources, 'sourceQuantities', ModuleName )
    call Allocate_test ( precisionVectors, noSources, 'precisionVectors', ModuleName )
    call Allocate_test ( precisionQuantities, noSources, 'precisionQuantities', ModuleName )
    ! Go round again and identify each quantity, work out what kind of file
    ! we're talking about
    precisionVectors = 0
    precisionQuantities = 0
    source = 0
    do keyNo = 2, nsons(node)
      l2gp_Version = 1
      son = subtree ( keyNo, node )
      fieldIndex = get_field_id ( son )
      select case ( fieldIndex )
      case ( f_source )
        source = source + 1
        gson = subtree(2,son)
        sourceVectors(source) = decoration(decoration(subtree(1,gson)))
        sourceQuantities(source) = decoration(decoration(decoration(subtree(2,gson))))
      case ( f_precision )
        if ( all ( outputType /= (/ l_l2gp, l_l2dgg /) ) ) &
          & call Announce_Error ( son, no_error_code, &
          & "Precision only appropriate for l2gp files" )
        gson = subtree(2,son)
        precisionVectors(source) = decoration(decoration(subtree(1,gson)))
        precisionQuantities(source) = decoration(decoration(decoration(subtree(2,gson))))
      case default
      end select
    end do

    if ( .not. checkpaths ) then
    ! Now go through and do some sanity checking
    do source = 1, noSources
      qty => GetVectorQtyByTemplateIndex ( vectors(sourceVectors(source)), &
      & sourceQuantities(source) )
      if ( qty%label == 0 ) call Announce_Error ( son, no_error_code, &
      & "Quantity does not have a label" )
      if ( precisionVectors(source) /= 0 ) then
      precQty => GetVectorQtyByTemplateIndex ( vectors(precisionVectors(source)), &
        & precisionQuantities(source) )
      ! Check that this is compatible with it's value quantitiy
      if ( qty%template%name /= precQty%template%name ) &
        & call Announce_Error ( son, no_error_code, &
        & "Precision and quantity do not match" )
      else
        precQty => NULL()
      end if
      ! Now check that things make sense
      if ( ValidateVectorQuantity ( qty, &
        & coherent=.true., stacked=.true., regular=.true., &
        & verticalCoordinate = (/ l_pressure, l_zeta, l_none/), &
        & minorFrame=.false., majorFrame=.false. ) ) then
        expectedType = l_l2gp
      else
        expectedType = l_l2aux
      end if
      if ( outputType /= expectedType .and. .not. &
        & ( outputType == l_l2dgg .and. expectedType == l_l2gp ) &
        &                           .and. .not.  &
        & ( outputType == l_l2fwm .and. expectedType == l_l2aux ) ) then
        call output ( "Offending quantity " )
        call display_string ( qty%template%name, strip=.true., advance='yes' )
        call Announce_Error ( son, no_error_code, &
          & "Inappropriate quantity for this file type in direct write" )
      end if
    end do
    end if
    
    ! Bail out at this stage if there is some kind of error.
    if ( error /= 0 ) return

    if ( DeeBUG ) then
      call output('Direct write to file', advance='yes')
      call output('File name: ', advance='no')
      call output(trim(filename), advance='yes')
      call output('hdfVersion: ', advance='no')
      call output(hdfVersion, advance='yes')
      call output('Num sources: ', advance='no')
      call output(noSources, advance='yes')
    endif

    ! If this is the first pass through, then we just log our request
    ! with the master
    if ( myMakeRequest ) then
      call LogDirectWriteRequest ( file, node )
    else
      ! OK, it's time to write this bit of the file
      if ( parallel%slave ) then
        createFileFlag = .false.
        if ( present ( create ) ) createFileFlag = create
      else
        createFileFlag = .not. any ( createdFilenames == file )
        if ( createFileFlag ) then
          noCreatedFiles = noCreatedFiles + 1
          if ( noCreatedFiles > maxFiles ) call MLSMessage ( &
            & MLSMSG_Error, ModuleName, 'Too many direct write files' )
          createdFilenames ( noCreatedFiles ) = file
        end if
      end if
      
      ! Open/create the file of interest
      call split_path_name(filename, path, file_base)
      if ( .not. TOOLKIT ) then
        handle = 0
        returnStatus = 0
      elseif ( any ( outputType == (/ l_l2gp /) ) ) then
        Handle = GetPCFromRef(file_base, mlspcf_l2gp_start, &
          & mlspcf_l2gp_end, &
          & TOOLKIT, returnStatus, l2gp_Version, DEEBUG, &
          & exactName=Filename)
      elseif ( any ( outputType == (/ l_l2dgg /) ) ) then
        Handle = GetPCFromRef(file_base, mlspcf_l2dgg_start, &
          & mlspcf_l2dgg_end, &
          & TOOLKIT, returnStatus, l2gp_Version, DEEBUG, &
          & exactName=Filename)
      elseif ( any ( outputType == (/ l_l2fwm /) ) ) then
        Handle = GetPCFromRef(file_base, mlspcf_l2fwm_full_start, &
          & mlspcf_l2fwm_full_end, &
          & TOOLKIT, returnStatus, l2gp_Version, DEEBUG, &
          & exactName=Filename)
      else
        Handle = GetPCFromRef(file_base, mlspcf_l2dgm_start, &
          & mlspcf_l2dgm_end, &
          & TOOLKIT, returnStatus, l2gp_Version, DEEBUG, &
          & exactName=Filename)
      end if
      if ( returnStatus /= 0 ) call MLSMessage ( &
         & MLSMSG_Error, ModuleName, &
         & 'Failed in GetPCFromRef for ' // trim(filename) )

      ! Done what we wished to do if just checking paths
      if ( checkPaths ) return
      
      if ( createFileFlag ) then
        fileaccess = DFACC_CREATE
      else
        fileaccess = DFACC_RDWR
      endif
      select case ( outputType )

      case ( l_l2gp, l_l2dgg )
        ! Before opening file, see which swaths are already there
        ! and which ones need to be created
        if ( DeeBUG ) print *, 'Allocating ', noSources
        nullify(createThisSource, nameBuffer)
        call Allocate_test ( createThisSource, noSources, 'createThisSource', &
          & ModuleName )
        call Allocate_test ( nameBuffer, noSources, 'nameBuffer', &
          & ModuleName )
        if ( .not. createFileFlag ) then
          do source = 1, noSources
            if(DEEBUG)print*,'Source:', source
            qty => GetVectorQtyByTemplateIndex ( vectors(sourceVectors(source)), &
              & sourceQuantities(source) )
            hdfNameIndex = qty%label
            if(DEEBUG)call display_string ( hdfNameIndex, strip=.true., advance='yes' )
            call get_string ( hdfNameIndex, nameBuffer(source), strip=.true. )
            if(DEEBUG)print*,'Done'
          enddo
          dummy = MLS_SWATH_IN_FILE(trim(fileName), nameBuffer, HdfVersion, &
            & createThisSource )
          if(DEEBUG)print*,'Got out of MLS_SWATH_IN_FILE'
        else
          createThisSource = .false.
        end if
        ! Call the l2gp open/create routine.  Filename is 'filename'
        ! file id should go into 'handle'
        handle = mls_io_gen_openF('sw', .true., ErrorType, &
          & record_length, FileAccess, trim(FileName), hdfVersion=hdfVersion)
      case ( l_l2aux, l_l2fwm )
        ! Call the l2aux open/create routine.  Filename is 'filename'
        ! file id should go into 'handle'
        handle = mls_io_gen_openF('hg', .true., ErrorType, &
          & record_length, FileAccess, FileName, hdfVersion=hdfVersion)
      end select
      
      if ( ErrorType /= 0 ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'DirectWriteCommand unable to open ' // trim(filename) )
      endif
      call SetUpNewDirect(newDirect, noSources)
      ! Add it to the database of directly writeable quantities
      dbindex = AddDirectToDatabase ( DirectDatabase, newDirect )
      call decorate ( node, -dbindex ) ! So we can find it later
      thisDirect => DirectDatabase(dbindex)
      ! Loop over the quantities to output
      do source = 1, noSources
        qty => GetVectorQtyByTemplateIndex ( vectors(sourceVectors(source)), &
          & sourceQuantities(source) )
        hdfNameIndex = qty%label
        call get_string ( hdfNameIndex, hdfName, strip=.true. )
        thisDirect%sdNames(source) = hdfName
        if ( precisionVectors(source) /= 0 ) then
          precQty => GetVectorQtyByTemplateIndex ( vectors(precisionVectors(source)), &
            & precisionQuantities(source) )
          ! Check that this is compatible with it's value quantitiy
          if ( qty%template%name /= precQty%template%name ) &
            & call Announce_Error ( son, no_error_code, &
            & "Precision and quantity do not match" )
        else
          precQty => NULL()
        end if
        
        if ( DeeBUG ) then
          call output('CreateFileFlag: ', advance='no')
          call output(createFileFlag, advance='yes')
          call output('file access: ', advance='no')
          call output(fileaccess, advance='yes')
          call output('file handle: ', advance='no')
          call output(handle, advance='yes')
          call output('outputType: ', advance='no')
          call output(outputType, advance='yes')
          call output('sd name: ', advance='no')
          call output(trim(hdfname), advance='yes')
        endif
        
        select case ( outputType )
        case ( l_l2gp, l_l2dgg )
          ! Call the l2gp swath write routine.  This should write the 
          ! non-overlapped portion of qty (with possibly precision in precQty)
          ! into the l2gp swath named 'hdfName' starting at profile 
          ! qty%template%instanceOffset + 1
          ! May optionally supply first, last profiles
          call DirectWrite_l2GP ( handle, qty, precQty, hdfName, chunkNo, &
            & hdfVersion, filename=filename, &
            & createSwath=(.not. createThisSource(source)) )
          if ( outputType == l_l2dgg ) then
            filetype=l_l2dgg
          else
            filetype=l_swath
          endif
        case ( l_l2aux, l_l2fwm )
          ! Call the l2aux sd write routine.  This should write the 
          ! non-overlapped portion of qty (with possibly precision in precQty)
          ! into the l2aux sd named 'hdfName' starting at profile 
          ! qty%template%instanceOffset ( + 1 ? )
          ! Note sure about the +1 in this case, probably depends whether it's a
          ! minor frame quantity or not.  This mixed zero/one indexing is beomming
          ! a real pain.  I wish I never want down that road!
          call DirectWrite_L2Aux ( handle, qty, precQty, hdfName, hdfVersion, &
            & chunkNo, chunks )
          filetype=l_hdf
        end select
      end do ! End loop over swaths/sds

      thisDirect%type = outputType
      thisDirect%fileNameBase = file_base

      ! Close the output file of interest (does this need to be split like this?)
      select case ( outputType )
      case ( l_l2gp, l_l2dgg )
        if ( DeeBUG ) print *, 'Deallocating ', noSources
        call Deallocate_test ( createThisSource, 'createThisSource', ModuleName )
        call Deallocate_test ( nameBuffer, 'nameBuffer', ModuleName )
        ! Call the l2gp close routine
        errortype = mls_io_gen_closeF('sw', Handle, hdfVersion=hdfVersion)
        ! errortype = he5_SWclose(Handle)
        if ( DeeBUG ) then
          print *, 'Tried to close ', trim(FIleName)
          print *, 'Handle ', Handle
          print *, 'hdfVersion ', hdfVersion
          print *, 'errortype ', errortype
        endif
        if ( hdfVersion == HDFVERSION_5 .and. FORCEDIRWRITEREOPEN ) then
          print *, 'Now forcibly opening and closing the hdf5'
          Handle =   mls_sfstart(trim(FileName), DFACC_RDONLY, hdfVersion, &
           & addingmetadata=.false. )
          if ( HANDLE == -1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'DirectWriteCommand unable to h5fopen ' // trim(filename) )
          errortype =   mls_sfend(Handle, hdfVersion, &
           & addingmetadata=.false. )
          if ( errortype /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'DirectWriteCommand unable to h5fclose ' // trim(filename) )
        endif
      case ( l_l2aux, l_l2fwm )
        ! Call the l2aux close routine
        errortype = mls_io_gen_closeF('hg', Handle, hdfVersion=hdfVersion)
      end select
      if ( errortype /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'DirectWriteCommand unable to close ' // trim(filename) )
      if ( createFileFlag .and. TOOLKIT .and. .not. SKIPMETADATA .and. &
        & outputType /= l_l2fwm ) then
        call add_metadata ( file_base, noSources, thisDirect%sdNames, &
          & hdfVersion, filetype, errortype )
        if ( errortype /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'DirectWriteCommand unable to addmetadata to ' // trim(filename) )
      endif
      
      ! Tell the master we're done
      if ( parallel%slave ) call FinishedDirectWrite ( ticket )
    end if

    call Deallocate_test ( sourceVectors, 'sourceVectors', ModuleName )
    call Deallocate_test ( sourceQuantities, 'sourceQuantities', ModuleName )
    call Deallocate_test ( precisionVectors, 'precisionVectors', ModuleName )
    call Deallocate_test ( precisionQuantities, 'precisionQuantities', ModuleName )

  end subroutine DirectWriteCommand

  ! ------------------------------------------------ LabelVectorQuantity -----
  subroutine LabelVectorQuantity ( node, vectors )
    use VectorsModule, only: VECTOR_T, VECTORVALUE_T, GETVECTORQTYBYTEMPLATEINDEX
    use MoreTree, only: GET_FIELD_ID, GET_BOOLEAN
    use Init_tables_module, only: F_QUANTITY, F_PREFIXSIGNAL, F_LABEL
    use Symbol_Table, only: ENTER_TERMINAL
    use Symbol_Types, only: T_STRING
    use String_Table, only: GET_STRING
    use MLSSignals_m, only: GETSIGNALNAME
    use Tree, only: NSONS, SUBTREE, SUB_ROSA, DECORATION
    ! Dummy arguments
    integer, intent(in) :: NODE          ! Tree node for l2cf line
    type (Vector_T), dimension(:), pointer :: VECTORS ! Vectors database
    ! Local variables
    integer :: FIELDINDEX               ! Type of field
    integer :: KEYNO                    ! Field index
    integer :: LABEL                    ! String index
    integer :: QUANTITYINDEX            ! Index into quantities database
    integer :: SON                      ! Tree node
    integer :: SOURCE                   ! Tree node
    integer :: VECTORINDEX              ! Index into database
    logical :: PREFIXSIGNAL             ! From l2cf
    ! Executable code
    type (VectorValue_T), pointer :: QTY ! The quantity
    character(len=1024) :: LABELSTR     ! The label itself

    prefixSignal = .false.
    ! Loop over the fields of the mlscf line
    do keyNo = 2, nsons(node) ! Skip spec name
      son = subtree(keyNo,node)
      fieldIndex = get_field_id(son)
      select case ( fieldIndex )
      case ( f_quantity )
        source = subtree(2,son) ! required to be an n_dot vertex
        vectorIndex = decoration(decoration(subtree(1,source)))
        quantityIndex = decoration(decoration(decoration(subtree(2,source))))
      case ( f_prefixSignal )
        prefixSignal = get_boolean ( son )
      case ( f_label )
        label = sub_rosa(subtree(2,son))
      case default ! Can't get here if tree_checker worked properly
      end select
    end do

    ! Get the quantity
    qty => GetVectorQtyByTemplateIndex ( vectors(vectorIndex), quantityIndex )

    ! Adapt the label if the prefix signal flag is set.
    if ( prefixSignal ) then
      if ( qty%template%signal == 0 ) then
        call Announce_Error ( node, no_error_code, &
          & 'The quantity has no signal so prefixSignal is not appropriate' )
        return
      end if
      call GetSignalName ( qty%template%signal, labelStr, &
        & sideband=qty%template%sideband )
      call Get_String( label, labelStr(len_trim(labelStr)+1:), strip=.true. )
      ! Now get an index for this possibly new name which may include the signal
      label = enter_terminal ( trim(labelStr), t_string, caseSensitive=.true. )
    end if

    ! Attach the label
    qty%label = label
  end subroutine LabelVectorQuantity

  ! --------------------------------------------------  JoinQuantities  -----
  ! This routine parses a line of the l2cf that is designed to join
  ! quantities together into l2gp/l2aux files
  subroutine JoinQuantities ( node, vectors, l2gpDatabase, l2auxDatabase, &
    & chunkNo, chunks )

    use Expr_m, only: EXPR
    use INIT_TABLES_MODULE, only: &
      & F_COMPAREOVERLAPS, F_FILE, F_HDFVERSION, F_OUTPUTOVERLAPS, &
      & F_PRECISION, F_PREFIXSIGNAL, F_SOURCE, F_SDNAME, F_SWATH, FIELD_FIRST, &
      & FIELD_LAST
    use INIT_TABLES_MODULE, only: L_PRESSURE, &
      & L_TRUE, L_ZETA, S_DIRECTWRITE, S_L2AUX, S_L2GP, S_TIME, S_LABEL
    use intrinsic, only: L_NONE, L_GEODANGLE, &
      & L_MAF, PHYQ_DIMENSIONLESS
    use L2AUXData, only: L2AUXData_T
    use L2GPData, only: L2GPData_T
    use L2ParInfo, only: PARALLEL, SLAVEJOIN
    use MLSCommon, only: MLSCHUNK_T, R8
    use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_M, only: GetSignalName
    use MoreTree, only: GET_BOOLEAN, GET_FIELD_ID, GET_SPEC_ID
    use String_Table, only: GET_STRING
    use Symbol_Table, only: ENTER_TERMINAL
    use Symbol_Types, only: T_STRING
    use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, NULL_TREE, SOURCE_REF, &
      & SUB_ROSA, SUBTREE
    use TREE_TYPES, only: N_NAMED, N_SET_ONE
    use VectorsModule, only: GetVectorQtyByTemplateIndex, &
      & ValidateVectorQuantity, Vector_T, VectorValue_T

    ! Dummy arguments
    integer, intent(in) :: NODE         ! The start of the l2cf line
    type (Vector_T), dimension(:), pointer :: vectors
    type (L2GPData_T), dimension(:), pointer :: l2gpDatabase
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase
    integer, intent(in) :: chunkNo
    type (MLSChunk_T), dimension(:), intent(in) :: chunks

    ! Local variables
    logical :: COMPAREOVERLAPS
    integer :: FIELDINDEX              ! F_..., see Init_Tables_Module
    integer :: FILE                 ! Name of output file for direct write
    integer :: SON                      ! Son of Key
    integer :: HDFVERSION               ! Version of hdf for directwrite
    integer :: HDFNAMEINDEX             ! Name of swath/sd
    integer :: KEY                      ! Index of an L2GP or L2AUX tree
    integer :: KEYNO                    ! Index of subtree of KEY
    integer :: MLSCFLine
    logical :: OutputOverlaps
    integer :: NAME                     ! Sub-rosa index of name of L2GP or L2AUX
    integer :: SOURCE                   ! Index in AST
    integer :: VALUE                    ! Value of a field
    integer :: VECTORINDEX              ! Index for vector to join
    integer :: QUANTITYINDEX            ! ind in qty tmpl database, not vector
    logical :: PREFIXSIGNAL             ! Prefix (i.e. make) the sd name the signal
    integer :: PRECVECTORINDEX          ! Index for precision vector
    integer :: PRECQTYINDEX             ! Index for precision qty (in database not vector)
    logical :: TIMING
    integer :: EXPRUNITS(2)                 ! From expr
    real (r8) :: EXPRVALUE(2)               ! From expr

    real :: T1, T2     ! for timing

    character(len=132) :: HDFNAME          ! Name for swath/sd
    logical :: GOT(field_first:field_last)
    type (VectorValue_T), pointer :: Quantity
    type (VectorValue_T), pointer :: PrecisionQuantity

    ! We know this node is named
    key = subtree(2,node)
    name = sub_rosa(subtree(1,node))

    got = .false.
    source = null_tree
    compareOverlaps = .false.
    outputOverlaps = .false.
    hdfNameIndex=name
    prefixSignal = .false.
    hdfVersion = 4

    ! Loop over the fields of the mlscf line
    do keyNo = 2, nsons(key) ! Skip spec name
      son = subtree(keyNo,key)
      fieldIndex = get_field_id(son)
      got(fieldIndex) = .true.
      select case ( fieldIndex )
      case ( f_source )
        source = subtree(2,son) ! required to be an n_dot vertex
        vectorIndex = decoration(decoration(subtree(1,source)))
        quantityIndex = decoration(decoration(decoration(subtree(2,source))))
      case ( f_precision )
        source = subtree(2,son) ! required to be an n_dot vertex
        precVectorIndex = decoration(decoration(subtree(1,source)))
        precQtyIndex = decoration(decoration(decoration(subtree(2,source))))
      case ( f_hdfVersion )
        call expr ( subtree(2,son), exprUnits, exprValue )
        if ( exprUnits(1) /= phyq_dimensionless ) &
          & call Announce_error ( son, NO_ERROR_CODE, &
          & 'No units allowed for hdfVersion: just integer 4 or 5')
        hdfVersion = exprValue(1)
      case ( f_prefixSignal )
        prefixSignal = get_boolean(son)
      case ( f_compareoverlaps )
        compareOverlaps = get_boolean(son)
      case ( f_outputoverlaps )
        outputOverlaps = get_boolean(son)
      case ( f_swath )
        hdfNameIndex = sub_rosa(subtree(2,son))
      case ( f_sdName )
        hdfNameIndex = sub_rosa(subtree(2,son))
      case ( f_file )
        file = sub_rosa(subtree(2,son))
      case default ! Can't get here if tree_checker worked properly
      end select
    end do
    
      ! Some final checks
    if ( any ( got ( (/ f_file, f_hdfVersion /) ) ) ) &
      & call Announce_Error ( key, NO_ERROR_CODE, &
      & 'File or hdfVersion not appropriate arguments for output l2aux/l2gp' )

    if ( error /= 0 ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "Errors in configuration prevent proceeding" )

    ! Identify the quantity
    quantity => GetVectorQtyByTemplateIndex(vectors(vectorIndex),quantityIndex)
    ! Get the precision quantity too perhaps
    if ( got ( f_precision ) ) then
      precisionQuantity => &
        & GetVectorQtyByTemplateIndex(vectors(precVectorIndex),precQtyIndex)
      if ( quantity%template%name /= precisionQuantity%template%name ) &
        & call announce_error(key, NO_ERROR_CODE, &
        & 'Quantity and precision quantity do not match')
    else
      precisionQuantity => NULL()
    end if
    
    ! Establish a swath/sd name for this quantity.
    hdfName = ''
    if ( prefixSignal ) &
      & call GetSignalName ( quantity%template%signal, hdfName, &
      &   sideband=quantity%template%sideband )
    call Get_String( hdfNameIndex, hdfName(len_trim(hdfName)+1:), strip=.true. )
    ! Now get an index for this possibly new name which may include the signal
    hdfNameIndex = enter_terminal ( trim(hdfName), t_string, caseSensitive=.true. )
    
    ! Now do the join, perhaps as a parallel slave, perhaps more directly.
    if ( parallel%slave ) then
      ! For slave tasks in a PVM system, simply ship this vector off to the master
      call SlaveJoin ( quantity, precisionQuantity, hdfName, key )
    else
      ! Now, depending on the properties of the source we deal with the
      ! vector quantity appropriately.
      if ( ValidateVectorQuantity ( quantity, &
        & coherent=.true., stacked=.true., regular=.true., &
        & minorFrame=.false., majorFrame=.false., &
        & verticalCoordinate = (/ l_pressure, l_zeta, l_none/) ) ) then 
        ! Coherent, stacked, regular quantities on pressure surfaces, or
        ! with no vertical coordinate system go in l2gp files.
        if ( get_spec_id(key) /= s_l2gp ) then
          call Announce_Error ( key, NO_ERROR_CODE, &
            & 'This quantity should be joined as an l2gp' )
          call MLSMessage ( MLSMSG_Error,&
            & ModuleName, 'This quantity should be joined as an l2gp')
        endif
        call JoinL2GPQuantities ( key, hdfNameIndex, quantity, &
          & precisionQuantity, l2gpDatabase, chunkNo )
      else
        ! All others go in l2aux files.
        if ( get_spec_id(key) /= s_l2aux ) then
          call Announce_Error ( key, NO_ERROR_CODE, &
            & 'This quantity should be joined as an l2aux' )
          call MLSMessage ( MLSMSG_Error,&
            & ModuleName, 'This quantity should be joined as an l2aux')
        endif
        call JoinL2AUXQuantities ( key, hdfNameIndex, quantity, &
          & l2auxDatabase, chunkNo, chunks )
      end if
    end if

  end subroutine JoinQuantities

  ! -----------------------------------------  JoinL2GPQuantities  -----

  ! This routine joins an l2gp line quantity to a database of such quantities.
  ! If this is the first time through, the database is created.

  ! The firstInstance and lastInstance arguments give an optional range of
  ! the instances that we wish to store in the l2gp quantity.  Otherwise, it
  ! defaults to the non overlapped region.

  subroutine JoinL2GPQuantities ( key, name, quantity, &
    & precision, l2gpDatabase, chunkNo, &
    & firstInstance, lastInstance, nameString )

    use INIT_TABLES_MODULE, only: L_PRESSURE, L_ZETA
    use intrinsic, only: L_NONE
    use L2GPData, only: AddL2GPToDatabase, ExpandL2GPDataInPlace, &
      & L2GPData_T, SetupNewL2GPRecord, RGP
    use MLSCommon, only: R4, R8, RV
    use String_Table, only: GET_STRING
    use TOGGLES, only: GEN, TOGGLE, LEVELS, SWITCHES
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, NULL_TREE, SOURCE_REF, &
      & SUB_ROSA, SUBTREE
    use VectorsModule, only: VectorValue_T

    ! Dummy arguments
    integer, intent(in) :: KEY          ! spec_args to Decorate with the L2GP index
    integer, intent(in) :: NAME         ! For the swath
    type (VectorValue_T), intent(in) :: QUANTITY ! Vector quantity
    type (VectorValue_T), pointer :: precision ! Optional vector quantity
    type (L2GPData_T), dimension(:), pointer :: L2GPDATABASE
    integer, intent(in) :: CHUNKNO
    integer, intent(in), optional :: FIRSTINSTANCE, LASTINSTANCE
    character(len=*), intent(in), optional :: nameString
    ! The last two are set if only part (e.g. overlap regions) of the quantity
    ! is to be stored in the l2gp data.

    ! Local variables
    type (L2GPData_T) :: NewL2GP
    type (L2GPData_T), pointer :: ThisL2GP
    integer :: Index
    integer :: FirstProfile, LastProfile ! Profile range in the l2gp to output to
    integer :: NoSurfsInL2GP, NoFreqsInL2GP
    integer :: UseFirstInstance, UseLastInstance, NoOutputInstances
    logical :: L2gpDataIsNew
    real(rv) :: HUGERGP
!   real(r8), dimension(:,:), pointer :: Values !??? Not used ???
    
    hugeRgp = real ( huge(0.0_rgp), rv )
    if ( toggle(gen) .and. levels(gen) > 0 ) &
      & call trace_begin ( "JoinL2GPQuantities", key )

    ! If this is the first chunk, we have to setup the l2gp quantity from
    ! scratch.  Otherwise, we expand it and fill up our part of it.
    l2gpDataIsNew = (.not. associated(l2gpDatabase))
    if ( .not. l2gpDataIsNew ) then
      index = decoration(key)
      l2gpDataIsNew = (index>=0)
    end if

    ! Work out what to do with the first and last instance information
    
    if ( present(firstInstance) ) then
      useFirstInstance = firstInstance
    else
      useFirstInstance = quantity%template%noInstancesLowerOverlap+1
    end if

    if ( present(lastInstance) ) then
      useLastInstance = lastInstance
    else
      useLastInstance = quantity%template%noInstances- &
        & quantity%template%noInstancesUpperOverlap
    end if

    noOutputInstances = useLastInstance-useFirstInstance+1
    ! If we've not been asked to output anything then don't carry on
    if ( noOutputInstances < 1 ) return

    if ( l2gpDataIsNew ) then
      ! Now create an empty L2GP record with this dimension

      if (any(quantity%template%verticalCoordinate == (/l_Pressure, l_Zeta /) )) then
        noSurfsInL2GP = quantity%template%noSurfs
      else
        noSurfsInL2GP = 0
      end if

      if ( quantity%template%frequencyCoordinate == l_None) then
         noFreqsInL2GP=0
      else
         noFreqsInL2GP=quantity%template%noChans
      end if

      call SetupNewL2GPRecord ( newL2GP, noFreqsInL2GP, noSurfsInL2GP )
      ! Setup the standard stuff, only pressure as it turns out.
      if ( quantity%template%verticalCoordinate == l_Pressure ) &
        & newL2GP%pressures = quantity%template%surfs(:,1)
      if ( quantity%template%verticalCoordinate == l_Zeta ) &
        & newL2GP%pressures = 10.0**(-quantity%template%surfs(:,1))
      ! It inherits its quantity type from the quantity template
      newL2GP%quantityType=quantity%template%quantityType
      ! Do something about frequency
      if ( associated ( quantity%template%frequencies ) ) then
        newL2GP%frequency = quantity%template%frequencies
      else
        newL2GP%frequency = 0.0
      end if

      ! Add it to the database of l2gp quantities
      index = AddL2GPToDatabase ( l2gpDatabase, newL2GP )
      ! Setup the pointer and index to be used later
      call decorate ( key, -index ) ! Remember where it is
      thisL2GP => l2gpDatabase(index)

    else
      ! Setup the index and pointer
      thisL2GP => l2gpDatabase(-index)
    end if

    ! Expand l2gp (initially all zero-size arrays) to take the new information
    call ExpandL2GPDataInPlace ( thisL2GP, &
      &thisL2GP%nTimes+noOutputInstances )

    ! Now copy the information from the quantity to the l2gpData

    ! name is an integer, but L2GP%name is Character data
    thisL2GP%nameIndex = name
    if ( present(nameString) ) then
      thisL2GP%name = nameString
    else
      call Get_String( name, thisL2GP%name, strip=.true.)
    end if
    lastProfile=thisL2GP%nTimes
    firstProfile=lastProfile-noOutputInstances+1

    ! Now fill the data, first the geolocation
    thisL2GP%latitude(firstProfile:lastProfile) = &
      & quantity%template%geodLat(1,useFirstInstance:useLastInstance)
    thisL2GP%longitude(firstProfile:lastProfile) = &
      & quantity%template%lon(1,useFirstInstance:useLastInstance)
    thisL2GP%solarTime(firstProfile:lastProfile) = &
      & quantity%template%solarTime(1,useFirstInstance:useLastInstance)
    thisL2GP%solarZenith(firstProfile:lastProfile) = &
      & quantity%template%solarZenith(1,useFirstInstance:useLastInstance)
    thisL2GP%losAngle(firstProfile:lastProfile) = &
      & quantity%template%losAngle(1,useFirstInstance:useLastInstance)
    thisL2GP%geodAngle(firstProfile:lastProfile) = &
      & quantity%template%phi(1,useFirstInstance:useLastInstance)
    thisL2GP%time(firstProfile:lastProfile) = &
      & quantity%template%time(1,useFirstInstance:useLastInstance)
    thisL2GP%chunkNumber(firstProfile:lastProfile)=chunkNo

    ! Now the various data quantities.

    ! For v0.1 we're only thinking about value.  The precision will
    ! come from matrices later in 0.5, and the diagnostics such as status
    ! and quality will come later too (probably 0.5, but maybe 1.0)

    thisL2GP%l2gpValue(:,:,firstProfile:lastProfile) = &
      & reshape ( max ( -hugeRgp, min ( hugeRgp, &
      &   quantity%values(:,useFirstInstance:useLastInstance) ) ), &
      &  (/max(thisL2GP%nFreqs,1),max(thisL2GP%nLevels,1),lastProfile-firstProfile+1/))
    if (associated(precision)) then
      thisL2GP%l2gpPrecision(:,:,firstProfile:lastProfile) = &
        & reshape ( max ( -hugeRgp, min ( hugeRgp, &
        &   precision%values(:,useFirstInstance:useLastInstance) ) ), &
        &  (/max(thisL2GP%nFreqs,1),max(thisL2GP%nLevels,1),lastProfile-firstProfile+1/))
    else
      thisL2GP%l2gpPrecision(:,:,firstProfile:lastProfile) = 0.0
    end if
    thisL2GP%status(firstProfile:lastProfile)='G'
    thisL2GP%quality(firstProfile:lastProfile)=0.0

    if ( toggle(gen) .and. levels(gen) > 0 ) call trace_end ( "JoinL2GPQuantities" )
  end subroutine JoinL2GPQuantities

  ! ----------------------------------------  JoinL2AUXQuantities  -----

  ! This subroutine is like the one above, except that the quantities it joins
  ! are destined to go in L2AUX quantities.

  subroutine JoinL2AUXQuantities ( key, name, quantity, l2auxDatabase, &
   & chunkNo, chunks, firstInstance, lastInstance )

    use intrinsic, only: L_NONE, L_GEODANGLE, &
      & L_MAF, PHYQ_DIMENSIONLESS
    use L2AUXData, only: AddL2AUXToDatabase, ExpandL2AUXDataInPlace, &
      & L2AUXData_T, L2AUXRank, SetupNewL2AUXRecord
    use MLSCommon, only: MLSCHUNK_T, R4, R8, RV
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use OUTPUT_M, only: BLANKS, OUTPUT
    use String_Table, only: DISPLAY_STRING, GET_STRING
    use TOGGLES, only: GEN, TOGGLE, LEVELS, SWITCHES
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, NULL_TREE, SOURCE_REF, &
      & SUB_ROSA, SUBTREE
    use VectorsModule, only: VectorValue_T

    ! Dummy arguments
    integer, intent(in) :: KEY     ! spec_args to decorate with the L2AUX index
    integer, intent(in) :: NAME    ! for the sd
    type (VectorValue_T), intent(in) :: quantity
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase
    integer, intent(in) :: chunkNo
    type (MLSChunk_T), dimension(:), intent(in) :: chunks
    integer, intent(in), optional :: firstInstance, lastInstance
    ! The last two are set if only part (e.g. overlap regions) of the quantity
    ! is to be stored in the l2aux data.

    ! Local variables
    logical :: DEEBUG
    integer :: FIRSTMAF
    integer :: FIRSTPROFILE
    integer :: DB_INDEX
    integer :: LASTMAF
    integer :: LASTPROFILE
    integer :: MAF
    integer :: NOMAFS
    integer :: NOOUTPUTINSTANCES
    integer :: USEFIRSTINSTANCE
    integer :: USELASTINSTANCE
    logical :: L2AUXDATAISNEW

    character (LEN=32) :: QUANTITYNAMESTR
    real(r8) :: HUGER4
    type (L2AUXData_T) :: NEWL2AUX
    type (L2AUXData_T), pointer :: THISL2AUX

    ! Executable code

    DEEBUG = (index(switches, 'join') /= 0)
    if ( toggle(gen) .and. levels(gen) > 0 ) &
      & call trace_begin ( "JoinL2AUXQuantities", key )

    hugeR4 = real ( huge(0.0_r4), r8 )

    if ( DEEBUG ) then
      call output('Joining vector quantity to L2AUX quantities', advance='yes')
      call output('minor frame? ', advance='no')
      call output(quantity%template%minorFrame, advance='no')
      call output('   major frame? ', advance='no')
      call output(quantity%template%majorFrame, advance='no')
      call output('   template  name ', advance='no')
      if ( quantity%template%name < 1 ) then
        call output('   (unnamed) ', advance='yes')
      else
        call display_string(quantity%template%name, strip=.true., advance='yes' )
      end if
    end if

    ! If this is the first chunk, we have to setup the l2aux quantity from
    ! scratch.  Otherwise, we expand it and fill up our part of it.
    l2auxDataIsNew = (.not. associated(l2auxDatabase))
    if ( .not. l2auxDataIsNew ) then
      db_index = decoration(key)
      l2auxDataIsNew = (db_index>=0)
    end if

    ! Work out what to do with the first and last Instance information
    if ( present(firstInstance) ) then
      useFirstInstance = firstInstance
    else
      useFirstInstance = quantity%template%noInstancesLowerOverlap+1
    end if
    if ( present(lastInstance) ) then
      useLastInstance = lastInstance
    else
      useLastInstance = quantity%template%noInstances- &
        & quantity%template%noInstancesUpperOverlap
    end if
    noOutputInstances = useLastInstance - useFirstInstance + 1

    ! If we've not been asked to output anything then don't carry on
    if ( noOutputInstances < 1 ) return

    if ( DEEBUG ) then
      call output('Joining L2Aux quantity with ', advance='no')
      call output(noOutputInstances, advance='no')
      call output(' instances ', advance='yes')
    end if

    ! Now if this is a new l2aux quantity, we need to setup an l2aux data type
    ! for it.
    if ( l2auxDataIsNew ) then
      ! We need to setup the quantity.  In some cases (minor/major frame)
      ! we can tell how big it is going to be.  Otherwise, create it empty
      if ( any ((/ quantity%template%minorFrame, quantity%template%majorFrame /)) ) then
        firstMAF = minval ( chunks%firstMAFIndex )
        lastMAF = maxval ( chunks%lastMAFIndex )
        noMAFs = lastMAF - firstMAF + 1
        if ( DEEBUG ) then
          call output('  firstMAF ', advance='no')
          call output(firstMAF, advance='no')
          call output('  noMAFs ', advance='no')
          call output(noMAFs, advance='yes')
        endif
      else
        ! Otherwise, we don't know how big it will be (at least in the Join
        ! scenario), so create it empty to begin with.
        firstMAF = 1
        noMAFs = 0
      end if
      ! Create the record accordingly
      call SetupNewL2AUXRecord ( newL2AUX, quantity%template, firstMAF, noMAFs )
      newL2AUX%instrumentModule=quantity%template%instrumentModule
      newL2AUX%quantityType=quantity%template%quantityType

      ! Add this l2aux to the database
      db_index = AddL2AUXToDatabase ( l2auxDatabase, newL2AUX )

      ! Setup the pointer and the index to be used later
      call decorate ( key, -db_index ) ! Remember where it is
      thisL2AUX => l2auxDatabase(db_index)
      thisL2AUX%name = name
    else
      ! Not a new l2aux, so just point ourselves to the old one.
      thisL2AUX => l2auxDatabase(-db_index)
    end if

    ! OK, now thisL2AUX points to an appropriate l2aux to fill, be it newly created,
    ! or old, and be it big enough or not.
    ! First, do we need to expand this?
    if ( .not. any ((/ quantity%template%minorFrame, quantity%template%majorFrame /)) ) then
      ! We need to expand this L2AUX to fit in the latest data.
      call ExpandL2AUXDataInPlace ( thisL2AUX, noOutputInstances )
    end if

    if ( quantity%template%minorFrame .or. quantity%template%majorFrame ) then
      ! Don't forget instanceOffset is for the first non-overlapped instance (ie MAF)
      ! Also remember the L2AUX data is already indexed from zero! Great!
      if (DEEBUG) call output ( "Doing the special calculation of first/last profile", advance='yes' )
      lastProfile = quantity%template%instanceOffset + quantity%template%noInstances - &
        & quantity%template%noInstancesUpperOverlap - &
        & quantity%template%noInstancesLowerOverlap - 1
    else
      lastProfile = thisL2AUX%dimensions(L2AUXRank)%noValues
    end if
    firstProfile = lastProfile - noOutputInstances + 1

    if ( DEEBUG ) then
      call output('  instance offset ' )
      call output( quantity%template%instanceOffset, advance='yes' )
      call output('  firstProfile ', advance='no')
      call output(firstProfile, advance='no')
      call output('   lastProfile ', advance='no')
      call output(lastProfile, advance='yes')
      call output('  FirstInstance ', advance='no')
      call output(useFirstInstance, advance='no')
      call output('   LastInstance ', advance='no')
      call output(useLastInstance, advance='yes')
      call output('  L2AUX%dimensions ', advance='no')
      call output(trim(thisL2AUX%dim_names), advance='yes')
      call output('  L2AUX%dim_units ', advance='no')
      call output(trim(thisL2AUX%dim_units), advance='yes')
      call output('  L2AUX%value_units ', advance='no')
      call output(trim(thisL2AUX%value_units), advance='yes')
      call output('  L2AUX%dimensions(1)%noValues ', advance='no')
      call output(thisL2AUX%dimensions(1)%noValues, advance='no')
      call output('  L2AUX%dimensions(2)%noValues ', advance='no')
      call output(thisL2AUX%dimensions(2)%noValues, advance='no')
      call output('  L2AUX%dimensions(3)%noValues ', advance='no')
      call output(thisL2AUX%dimensions(3)%noValues, advance='yes')
      call output('shape(l2aux values) ', advance='no')
      call output(shape(thisL2AUX%values), advance='yes')
      if ( any ( thisL2AUX%dimensions(L2AUXRank)%dimensionFamily &
        & == (/ L_GeodAngle, L_MAF /) ) ) then
        call output('   dimensions ', advance='no')
        call output(size(thisL2AUX%dimensions(L2AUXRank)%values), advance='no')
      else
        call output(' (dimensions unassociated)', advance='no')
      end if
      call output('   values 3rd coord ', advance='no')
      call output(size(thisL2AUX%values(1,1,:)), advance='yes')
    end if

    select case (thisL2AUX%dimensions(L2AUXRank)%dimensionFamily)
    case ( L_GeodAngle )
      thisL2AUX%dimensions(L2AUXRank)%values(firstProfile:lastProfile)=&
        & quantity%template%phi(1,useFirstInstance:useLastInstance)
    case ( L_MAF )
      do maf = firstProfile, lastProfile
        thisL2AUX%dimensions(L2AUXRank)%values(maf) = maf
      end do
    case default
    end select
    
    ! Check that reshape has a prayer of succeeding
    if ( DEEBUG ) then
      call output('  num l2aux values/profile ', advance='no')
      call output(size(thisL2AUX%values, 1)*size(thisL2AUX%values, 2), &
       & advance='yes')
      call output('   num dim values/profile ', advance='no')
      call output(thisL2AUX%dimensions(1)%noValues* &
       &              thisL2AUX%dimensions(2)%noValues, advance='yes')
      call output('  num l2aux values (total) ', advance='no')
      call output(size(thisL2AUX%values, 1)*size(thisL2AUX%values, 2)* &
       &              (lastProfile-firstProfile+1), advance='yes')
      call output('   num qty values ', advance='no')
      call output(size(quantity%values, 1)* &
       &              (useLastInstance-useFirstInstance+1), advance='yes')
    endif
    if ( size(thisL2AUX%values, 1)*size(thisL2AUX%values, 2) &
     & /= thisL2AUX%dimensions(1)%noValues*thisL2AUX%dimensions(2)%noValues ) &
     & call MLSMessage ( MLSMSG_Error, &
     & ModuleName, "Reshape fails: size mismatch betw. dims and values" )
    if ( size(thisL2AUX%values, 1)*size(thisL2AUX%values, 2)*(lastProfile-firstProfile+1) &
     & /= size(quantity%values, 1)*(useLastInstance-useFirstInstance+1) ) &
     & call MLSMessage ( MLSMSG_Error, &
     & ModuleName, "Reshape fails: size mismatch betw. quantity and l2aux" )
    thisL2AUX%values(:,:,firstProfile:lastProfile) = &
      & reshape ( max ( -hugeR4, min ( hugeR4, &
      & quantity%values(:,useFirstInstance:useLastInstance) ) ), &
      &   (/ thisL2AUX%dimensions(1)%noValues, &
      &      thisL2AUX%dimensions(2)%noValues, &
      &      lastProfile-firstProfile+1/) )
    
    if ( toggle(gen) .and. levels(gen) > 0 ) &
      & call trace_end ( "JoinL2AUXQuantities" )

  end subroutine JoinL2AUXQuantities

! =====     Private Procedures     =====================================

  ! ---------------------------------------------  Announce_Error  -----
  subroutine ANNOUNCE_ERROR ( where, CODE, ExtraMessage, FIELDINDEX )

    use intrinsic, only: FIELD_INDICES
    use LEXER_CORE, only: PRINT_SOURCE
    use OUTPUT_M, only: BLANKS, OUTPUT
    use String_Table, only: DISPLAY_STRING
    use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, NULL_TREE, SOURCE_REF, &
      & SUB_ROSA, SUBTREE

    integer, intent(in) :: where   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message
    integer, intent(in), optional :: FIELDINDEX ! Extra information for msg
    character (LEN=*), intent(in), optional :: ExtraMessage

    error = max(error,1)
    call output ( '***** At ' )
    if ( where > 0 ) then
      call print_source ( source_ref(where) )
    else
      call output ( '(no lcf tree available)' )
    end if
    call output ( ': ' )
    select case ( code )
      case ( NotAllowed )
        call output('Field ')
        call display_string(field_indices(fieldIndex))
        call output(' is not allowed in this context',advance='yes')
      case default
        call output ( " command caused an unrecognized programming error", advance='yes' )
    end select
    if ( present(ExtraMessage) ) then
      call output(ExtraMessage, advance='yes')
    end if
  end subroutine ANNOUNCE_ERROR

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Join

!
! $Log$
! Revision 2.96  2003/11/07 00:46:51  pwagner
! New quicker preflight option: --checkPaths
!
! Revision 2.95  2003/10/20 18:21:45  pwagner
! Timings breakdown added for directWrite
!
! Revision 2.94  2003/10/10 00:00:24  pwagner
! Should quit properly if SIPS and no PCFid match for file name
!
! Revision 2.93  2003/09/12 21:45:52  pwagner
! Only prints l2gp label during DirectWrite if DEEBUG
!
! Revision 2.92  2003/09/04 22:40:04  pwagner
! Gets dgg file name from dgg PCFid when DirectWrite usingPCF
!
! Revision 2.91  2003/09/03 23:05:49  livesey
! More problems with precision in DirectWrite, hope I've got them all now.
!
! Revision 2.90  2003/09/03 00:53:50  livesey
! Bug fix on directWrite, it was storing the value field in the precision
! too!
!
! Revision 2.89  2003/08/28 23:51:58  livesey
! Renamed some variables to make them more obvious
!
! Revision 2.88  2003/08/14 20:11:30  pwagner
! DirectWrite may take l2fwm types for fwm radiances
!
! Revision 2.87  2003/08/01 20:38:31  pwagner
! Distinguishes between l2dgg and l2gp when writing metadata as part of directwrite
!
! Revision 2.86  2003/07/25 00:51:06  livesey
! Added file type l2dgg to support metadata.
!
! Revision 2.85  2003/07/23 18:30:35  cvuu
! reduce routine printing
!
! Revision 2.84  2003/07/15 23:39:01  pwagner
! Disabled most printing
!
! Revision 2.83  2003/07/11 01:24:20  livesey
! More changes trying to get the direct write going.
!
! Revision 2.82  2003/07/09 21:49:53  pwagner
! Tries to figure out in advance whether to create swath
!
! Revision 2.81  2003/07/08 00:15:51  livesey
! Various tidy ups and reworks
!
! Revision 2.80  2003/07/07 23:52:13  pwagner
! Slave that creates DirectWrite file may also add_metadata
!
! Revision 2.79  2003/07/07 20:29:43  livesey
! Mainly cosmetic changes
!
! Revision 2.78  2003/07/07 17:31:11  livesey
! Various things to get DirectWrite working
!
! Revision 2.77  2003/07/02 00:55:27  pwagner
! Some improvements in DirectWrites of l2aux, l2gp
!
! Revision 2.76  2003/06/26 23:13:52  pwagner
! New debugging output; distinguishes between l2gp/l2aux quantities better
!
! Revision 2.75  2003/06/24 23:54:07  pwagner
! New db indexes stored for entire direct file
!
! Revision 2.74  2003/06/24 23:30:00  livesey
! Finished LabelVectorQuantity and made some other bug fixes.
!
! Revision 2.73  2003/06/23 23:55:17  pwagner
! Added DirectData_T to keep track of data written directly
!
! Revision 2.72  2003/06/20 19:38:25  pwagner
! Allows direct writing of output products
!
! Revision 2.71  2003/05/30 00:09:27  livesey
! Can now directWrite major frame quantities
!
! Revision 2.70  2003/05/12 02:06:23  livesey
! Bug fix for prefixSignal L2GPs and also bound r8->r4 conversion
!
! Revision 2.69  2003/02/08 00:31:31  pwagner
! Now saves quantityType in newl2gp
!
! Revision 2.68  2003/01/30 01:03:24  pwagner
! Stores quantity type taken from source vector in l2aux
!
! Revision 2.67  2003/01/17 23:11:26  pwagner
! Moved most ops out of LoinL2AUXData to SetupL2AUXData
!
! Revision 2.66  2002/12/19 15:53:47  livesey
! Allowed verticalCoordinate=l_none quantities back into the l2gp fold.
!
! Revision 2.65  2002/11/26 23:38:01  livesey
! Better joining of major frame quantities
!
! Revision 2.64  2002/10/29 21:54:21  livesey
! Made join less verbose.
!
! Revision 2.63  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.62  2002/08/20 22:10:50  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.61  2002/08/20 20:10:30  livesey
! Dealt with frequency in l2gps
!
! Revision 2.60  2002/05/28 17:09:57  livesey
! Removed print statements
!
! Revision 2.59  2002/05/22 00:48:52  livesey
! Added direct write stuff
!
! Revision 2.58  2002/05/16 22:36:46  livesey
! Fixed a bug with joining minor frame quantities with overlaps.
!
! Revision 2.57  2002/04/08 20:49:17  pwagner
! Swath name optionally passed to JoinL2GPQuantities
!
! Revision 2.56  2002/04/06 00:35:21  pwagner
! Should accept actual case of swathname for l2gp
!
! Revision 2.55  2002/03/20 00:46:47  pwagner
! Removed 2 unused lits
!
! Revision 2.54  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.53  2001/10/30 01:45:21  livesey
! Some modifications/fixes to parallel join
!
! Revision 2.52  2001/10/08 23:38:58  pwagner
! Tiny fixes; not perfect yet
!
! Revision 2.51  2001/10/06 00:27:42  pwagner
! Still some problems with diagnostics
!
! Revision 2.50  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.49  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
! Revision 2.48  2001/09/08 00:21:44  pwagner
! Revised to work for new column Abundance in lone swaths
!
! Revision 2.47  2001/09/05 20:34:56  pwagner
! Reverted to pre-columnAbundance state
!
! Revision 2.46  2001/08/03 23:13:52  pwagner
! Began testing; at least now exits normally again
!
! Revision 2.45  2001/08/02 23:58:31  pwagner
! More complete treatment of column abundance(s)
!
! Revision 2.44  2001/08/02 00:18:55  pwagner
! Began adding column quantities; incomplete
!
! Revision 2.43  2001/07/31 23:25:32  pwagner
! Able to accept 2 new fields for join of column; does nothing yet
!
! Revision 2.42  2001/06/19 22:52:31  pwagner
! l_none  no longer got from init_tables_module
!
! Revision 2.41  2001/05/23 21:59:43  livesey
! Interim version, almost there
!
! Revision 2.40  2001/05/23 01:43:19  livesey
! New parallel version in progress
!
! Revision 2.39  2001/05/19 01:19:58  livesey
! Fills precision correctly for l2gp!
!
! Revision 2.38  2001/05/14 22:23:53  livesey
! Embarassing bug fix, to do with renumbering of minor frame L2AUX quantities.
!
! Revision 2.37  2001/05/12 00:18:17  livesey
! Tidied up array bounds for L2AUX/MAF.
!
! Revision 2.36  2001/05/10 16:31:24  livesey
! Added prefix signal option for swath/sd name
!
! Revision 2.35  2001/05/08 23:25:32  livesey
! Added the precision stuff for l2gp's
!
! Revision 2.34  2001/05/08 21:51:02  livesey
! Removed some old xStar, yStar, kStar stuff.
!
! Revision 2.33  2001/05/03 20:32:19  vsnyder
! Cosmetic changes
!
! Revision 2.32  2001/05/02 22:22:43  pwagner
! Removed SDPToolkit use
!
! Revision 2.31  2001/04/28 01:30:14  livesey
! Basically gone back to an earlier version.  As l2pc's now output
! directly as matrices there is no need for Join to think about them.
!
! Revision 2.30  2001/04/27 21:52:39  livesey
! Removed l2pc stuff
!
! Revision 2.29  2001/04/26 20:02:09  livesey
! Made l2pc database a saved array in L2PC_m
!
! Revision 2.28  2001/04/26 15:59:30  livesey
! Tidied up uses
!
! Revision 2.27  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.26  2001/04/26 00:07:16  livesey
! Insulate vector is gone
!
! Revision 2.25  2001/04/25 21:54:22  livesey
! Added candol2pc flag
!
! Revision 2.24  2001/04/25 21:51:46  livesey
! Tidied up Join for l2pcs
!
! Revision 2.23  2001/04/25 20:33:28  livesey
! Minor improvements to Join l2pc stuff
!
! Revision 2.22  2001/04/24 20:20:27  livesey
! L2PC moved to lib and word bin dropped from types etc.
!
! Revision 2.21  2001/04/24 20:04:54  livesey
! Added l2pc joining
!
! Revision 2.20  2001/04/10 23:44:44  vsnyder
! Improve 'dump'
!
! Revision 2.19  2001/03/15 23:26:56  livesey
! Avoid calling ExpandL2AUXQuantitiesInPlace for minor frame quantities.
! Really saves on memory thrashing.
!
! Revision 2.18  2001/03/15 21:18:57  vsnyder
! Use Get_Spec_ID instead of decoration(subtree...
!
! Revision 2.17  2001/03/06 22:40:41  livesey
! New L2AUX stuff
!
! Revision 2.16  2001/03/05 20:46:41  livesey
! Removed a debugging statement left behind
!
! Revision 2.15  2001/03/05 01:19:45  livesey
! Removed a print statement
!
! Revision 2.14  2001/03/05 01:01:12  livesey
! Bug fix, now uses GetVectorQtyFromTemplateIndex
!
! Revision 2.13  2001/03/01 18:38:27  livesey
! Fixed bug with verticalCoordinate==l_Zeta
!
! Revision 2.12  2001/02/27 17:38:21  livesey
! Tidied things up, removed unnecessary arguments
!
! Revision 2.11  2001/02/27 00:50:31  livesey
! Added ability to Join verticalCoordinate=l_zeta quantities into l2gp entities.
!
! Revision 2.10  2001/02/16 00:50:17  livesey
! Added error to avoid confusion with L2GP in ReadApriori section
!
! Revision 2.9  2001/02/09 19:30:16  vsnyder
! Move checking for required and duplicate fields to init_tables_module
!
! Revision 2.8  2001/02/09 18:01:46  livesey
! Various further updates, set default values for status and quality
!
! Revision 2.7  2001/02/09 00:38:22  livesey
! Various updates
!
! Revision 2.6  2001/01/03 18:15:13  pwagner
! Changed types of t1, t2 to real
!
! Revision 2.5  2000/11/16 02:19:01  vsnyder
! Implement timing.
!
! Revision 2.4  2000/11/13 23:02:21  pwagner
! Adapted for rank2 vectorsModule
!
! Revision 2.3  2000/10/05 16:37:19  pwagner
! Now compiles with new L2GPData module
!
! Revision 2.2  2000/09/11 19:34:35  ahanzel
! Removed old log entries in file.
!
! Revision 2.1  2000/09/08 22:55:56  vsnyder
! Revised to use the tree output by the parser
!

