! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=======================================================================================

module OutputAndClose ! outputs all data from the Join module to the
                      ! appropriate L2 Files

!=======================================================================================

  use DirectWrite_m, only: DirectData_T
! use DirectWrite_m, only: dump
  use Hdf, only: DFACC_CREATE
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
  use OUTPUT_M, only: blanks, OUTPUT
  use STRING_TABLE, only: GET_STRING

  implicit none
  private
  public :: Output_Close, add_metadata

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = &
    & "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

  ! -----     Private declarations     ---------------------------------

  ! Should we get this from MLSL2Options?
  logical, parameter :: LOGFILEGETSMETADATA = .false.  ! metadata to log file?
  ! For Announce_Error
  integer :: ERROR

contains ! =====     Public Procedures     =============================

  ! -----------------------------------------------  Output_Close  -----
  subroutine Output_Close ( root, l2gpDatabase, l2auxDatabase, DirectDatabase, &
    & matrices, canWriteL2PC )

    ! Hard-wired assumptions:

    ! ----------------------- metadata ------------------------

    !   for the l2aux the mcf is mlspcf_mcf_l2dgm_start
    !   for the log file the mcf is mlspcf_mcf_l2log_start
    !   for the dgg file the mcf is mlspcf_mcf_l2dgg_start

    ! The correspondence between MCF and l2gp files is determined by
    ! the value of        MCFFORL2GPOPTION
    ! (see write_metadata module for fuller explanation)

    use Allocate_Deallocate, only: Deallocate_Test, Allocate_test
    use Expr_M, only: Expr
    use INIT_TABLES_MODULE, only: F_ASCII, F_FILE, F_HDFVERSION, &
      & F_METANAME, F_METADATAONLY, F_OVERLAPS, F_PACKED, F_QUANTITIES, &
      & F_TYPE, F_WRITECOUNTERMAF, F_DONTPACK, &
      & L_L2AUX, L_L2DGG, L_L2GP, L_L2PC, S_OUTPUT, S_TIME
    use Intrinsic, only: l_swath, l_hdf, PHYQ_Dimensionless
    use L2AUXData, only: L2AUXDATA_T, cpL2AUXData, WriteL2AUXData
    use L2GPData, only: AVOIDUNLIMITEDDIMS, L2GPNameLen, L2GPData_T, &
      & MAXSWATHNAMESBUFSIZE, cpL2GPData, WriteL2GPData
    use L2PC_m, only: WRITEONEL2PC, OUTPUTHDF5L2PC
    use MatrixModule_1, only: MATRIX_DATABASE_T, MATRIX_T, GETFROMMATRIXDATABASE
    use MLSCommon, only: I4
    use MLSFiles, only: HDFVERSION_5, &
      & GetPCFromRef, mls_exists, &
      & MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF, MLS_SFSTART, MLS_SFEND, &
      & SPLIT_PATH_NAME, unSplitName
    use MLSL2Options, only: CATENATESPLITS, CHECKPATHS, &
      & DEFAULT_HDFVERSION_WRITE, PENALTY_FOR_NO_METADATA, TOOLKIT
    use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES
    use MLSPCF2, only: MLSPCF_L2DGM_END, MLSPCF_L2DGM_START, MLSPCF_L2GP_END, &
      & MLSPCF_L2GP_START, mlspcf_l2dgg_start, mlspcf_l2dgg_end, &
      & Mlspcf_mcf_l2gp_start, Mlspcf_mcf_l2dgm_start, &
      & Mlspcf_mcf_l2dgg_start
    use MLSSets, only: FindFirst, FindNext
    use MLSStrings, only: Array2List
    use MoreTree, only: Get_Spec_ID, GET_BOOLEAN
    use SDPToolkit, only: PGS_S_SUCCESS, PGSD_IO_GEN_WSEQFRM, Pgs_smf_getMsg
    use Time_M, only: Time_Now
    use TOGGLES, only: GEN, TOGGLE, Switches
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATION, NODE_ID, NSONS, SUBTREE, SUB_ROSA
    use TREE_TYPES, only: N_NAMED
    use WriteMetadata, only: L2PCF, WriteMetaLog

    ! Arguments
    integer, intent(in) :: ROOT   ! Of the output section's AST
    type (L2GPData_T), dimension(:), pointer :: L2GPDATABASE ! L2GP products
    type (L2AUXData_T), dimension(:), pointer :: L2AUXDATABASE ! L2AUX products
    type (Matrix_Database_T), dimension(:), pointer :: MATRICES ! Matrix database (for l2pcs)
    type (DirectData_T), dimension(:), pointer :: DirectDatabase
    logical, intent(in) :: canWriteL2PC ! Flag
    ! type(PCFData_T), intent(in) :: l2pcf

    ! - - - Local declarations - - -

    logical :: ASCII                    ! Is this l2pc ascii?
    integer :: DB_index
    integer, dimension(:), pointer :: DONTPACK ! Quantities not to pack
    logical, parameter :: DEBUG = .FALSE.
    integer :: FIELD_INDEX              ! F_... field code
    integer :: FIELD_NO                 ! Index of assign vertex sons of Key
    integer :: FIELDVALUE               ! For get_boolean
    character (len=132) :: FILE_BASE    ! From the FILE= field
    integer :: GSON                     ! Son of Son -- an assign node
    integer :: hdfVersion               ! 4 or 5 (corresp. to hdf4 or hdf5)
    integer :: IN_FIELD_NO              ! Index of sons of assign vertex
    integer :: KEY                      ! Index of spec_args node
    integer :: L2auxFileHandle, L2aux_Version
    character (len=132) :: L2auxPhysicalFilename
    integer :: L2aux_mcf, L2dgg_mcf, L2gp_mcf  ! mcf numbers for writing metadata
    integer :: L2gpFileHandle, L2gp_Version
    character (len=132) :: L2gpPhysicalFilename
    integer :: L2PCUNIT
    logical :: madeFile
    integer, parameter:: MAXQUANTITIESPERFILE=10000
    integer :: Metadata_error
    character (len=32) :: meta_name    ! From the metaName= field
    character (len=32) :: Mnemonic
    character (len=256) :: Msg
    integer :: NAME                     ! string index of label on output
    integer :: NODE
    integer :: Numquantitiesperfile     ! < MAXQUANTITIESPERFILE
    integer :: OUTPUT_TYPE              ! L_L2AUX, L_L2GP, L_PC, L_L2DGG
    logical :: PACKED                   ! Do we pack this l2pc?
    character (len=132) :: path
    integer :: QUANTITIESNODE           ! A tree node
    character(len=L2GPNameLen), dimension(MAXQUANTITIESPERFILE) :: QuantityNames  ! From "quantities" field
    integer :: RECLEN                   ! For file stuff
    integer :: record_length
    integer :: ReturnStatus
    integer(i4) :: SDFID                ! File handle
    character (len=MAXSWATHNAMESBUFSIZE) :: sdList
    integer :: SON                      ! Of Root -- spec_args or named node
    integer :: SPEC_NO                  ! Index of son of Root
    integer :: SWFID
    real :: T1, T2     ! for timing
    integer :: Type                     ! Type of value returned by EXPR
    integer :: Units(2)                 ! Units of value returned by EXPR
    double precision :: Value(2)        ! Value returned by EXPR
    type (Matrix_T), pointer :: TMPMATRIX ! A pointer to a matrix to write into l2pc
    logical :: TIMING
    logical :: WriteCounterMAF          ! Add the counter MAF field
    logical :: WriteMetaDataOnly        ! Because it was a directWrite

    ! Executable code
    timing = section_times
    if ( timing ) call time_now ( t1 )
    nullify ( dontPack )

    if ( toggle(gen) ) call trace_begin ( "Output_Close", root)

    error = 0

    if (index(switches, 'pro') /= 0) then
      call output ( '============ Level 2 Products ============', advance='yes' )
      call output ( ' ', advance='yes' )
    end if

    ! l2gp_mcf will be incremented in get_l2gp_mcf (if MCFFORL2GPOPTION == 1)
    l2gp_mcf = mlspcf_mcf_l2gp_start - 1   

    l2aux_mcf = mlspcf_mcf_l2dgm_start
    l2dgg_mcf = mlspcf_mcf_l2dgg_start

    ! Loop over the lines in the l2cf

    do spec_no = 2, nsons(root)-1 ! Skip name at begin and end of section

      l2gp_Version = 1
      l2aux_Version = 1
      hdfVersion = DEFAULT_HDFVERSION_WRITE
      meta_name = ''
      writeCounterMAF = .false.
      writeMetaDataOnly = .false.

      son = subtree(spec_no,root)
      if ( node_id(son) == n_named ) then ! Is spec labeled?
        key = subtree(2,son)
        name = sub_rosa(subtree(1,son))
      else ! Son is n_spec_args
        key = son
        name = 0
      end if

      select case( get_spec_id(key) )
      case ( s_output )
        do field_no = 2, nsons(key)       ! Skip the command name
          gson = subtree(field_no, key)   ! An assign node
          if ( nsons(gson) > 1 ) then
            fieldValue = decoration(subtree(2,gson)) ! The field's value
          else
            fieldValue = gson
          end if
          field_index = decoration(subtree(1,gson))
          select case ( field_index )   ! Field name
          case ( f_file )
            call get_string ( sub_rosa(subtree(2,gson)), file_base )
            file_base = file_base(2:LEN_TRIM(file_base)-1) ! Parser includes quotes
          case ( f_metaName )
            call get_string ( sub_rosa(subtree(2,gson)), meta_name )
            meta_name = meta_name(2:LEN_TRIM(meta_name)-1) ! Parser includes quotes
          case ( f_type )
            output_type = decoration(subtree(2,gson))
          case ( f_writeCounterMAF )
            writeCounterMAF = get_boolean ( fieldValue )
          case ( f_MetaDataOnly )
            writeMetaDataOnly = get_boolean ( fieldValue )
          case ( f_hdfVersion )
            call expr ( subtree(2,gson), units, value, type )
            if ( units(1) /= phyq_dimensionless ) &
              & call Announce_error ( gson, &
              & 'No units allowed for hdfVersion: just integer 4 or 5')
            hdfVersion = value(1)
          case default                  ! Everything else processed later
          end select
        end do

        if ( DEBUG ) call output('l2gp type number: ', advance='no')
        if ( DEBUG ) call output(l_l2gp, advance='yes')

        if ( DEBUG ) call output('l2aux type number: ', advance='no')
        if ( DEBUG ) call output(l_l2aux, advance='yes')

        if ( DEBUG ) call output('l2dgg type number: ', advance='no')
        if ( DEBUG ) call output(l_l2dgg, advance='yes')

        if ( DEBUG ) call output('output type number: ', advance='no')
        if ( DEBUG ) call output(output_type, advance='yes')

        if ( DEBUG ) call output('file_base: ', advance='no')
        if ( DEBUG ) call output(trim(file_base), advance='yes')

        ! Otherwise--normal output commands
        select case ( output_type )
        case ( l_l2gp ) ! --------------------- Writing l2gp files -----
          if ( DEBUG ) call output('output file type l2gp', advance='yes')
          ! Get the l2gp file name from the PCF

          if ( TOOLKIT ) then
            call split_path_name(file_base, path, file_base)
            if ( DEBUG ) call output('file_base after split: ', advance='no')
            if ( DEBUG ) call output(trim(file_base), advance='yes')

            l2gpFileHandle = GetPCFromRef(file_base, mlspcf_l2gp_start, &
              & mlspcf_l2gp_end, &
              & TOOLKIT, returnStatus, l2gp_Version, DEBUG, &
              & exactName=l2gpPhysicalFilename)
          else
            l2gpPhysicalFilename = file_base
            returnStatus = 0
          end if

          if ( returnStatus == 0 .and. .not. checkPaths ) then
            if ( DEBUG ) call output(&
              & 'file name: ' // TRIM(l2gpPhysicalFilename), advance='yes')
            ! Open the HDF-EOS file and write swath data

            if ( DEBUG ) call output('Attempting swopen', advance='yes')
            !           swfid = swopen(l2gpPhysicalFilename, DFACC_CREATE)
            swfid = mls_io_gen_openF('swopen', .TRUE., returnStatus, &
              & record_length, DFACC_CREATE, FileName=l2gpPhysicalFilename, &
              & hdfVersion=hdfVersion, debugOption=.false. )

            ! Loop over the segments of the l2cf line

            numquantitiesperfile = 0
            do field_no = 2, nsons(key) ! Skip "output" name
              gson = subtree(field_no,key)
              select case ( decoration(subtree(1,gson)) )
              case ( f_quantities )
                do in_field_no = 2, nsons(gson)
                  db_index = -decoration(decoration(subtree(in_field_no ,gson)))
                  if ( db_index >= 1 ) then
                    call writeL2GPData ( l2gpDatabase(db_index), swfid, &
                      & hdfVersion=hdfVersion )
                    numquantitiesperfile = numquantitiesperfile+1
                    if ( numquantitiesperfile > MAXQUANTITIESPERFILE ) then
                      call announce_error ( son, &
                        & 'Attempt to write too many l2gp quantities to a file', &
                        & numquantitiesperfile )
                      numquantitiesperfile = MAXQUANTITIESPERFILE
                    end if
                    quantityNames(numquantitiesperfile) = l2gpDatabase(db_index)%name
                  else
                    call MLSMessage ( MLSMSG_Warning, ModuleName, &
                      & 'Unable to write quantity to l2gp file, perhaps no chunks processed' )
                  end if
                end do ! in_field_no = 2, nsons(gson)
              case ( f_overlaps )
                ! ??? More work needed here
              end select
            end do ! field_no = 2, nsons(key)

            if ( DEBUG ) call output('Attempting swclose', advance='yes')
            !           returnStatus = swclose(swfid)
            returnStatus = mls_io_gen_closeF('swclose', swfid, &
              & hdfVersion=hdfVersion)
            if ( returnStatus /= PGS_S_SUCCESS ) then
              call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
              call MLSMessage ( MLSMSG_Error, ModuleName, &
                &  "Error closing  l2gp file:  "//mnemonic//" "//msg )
            else if (index(switches, 'pro') /= 0) then
              call announce_success(l2gpPhysicalFilename, 'l2gp', &
                & numquantitiesperfile, quantityNames, hdfVersion=hdfVersion)
            end if

            if ( .not. TOOLKIT ) cycle

            ! Write the metadata file
            call add_metadata ( son, file_base, numquantitiesperfile, &
              & quantityNames, hdfVersion, l_swath, metadata_error )
          else if ( returnStatus /= PGS_S_SUCCESS ) then
            call announce_error ( son, &
              &  "Error finding l2gp file matching:  "//file_base, returnStatus)
          end if

        case ( l_l2aux ) ! ------------------------------ Writing l2aux files ---

          if ( DEBUG ) call output ( 'output file type l2aux', advance='yes' )
          ! Get the l2aux file name from the PCF

          if ( TOOLKIT ) then
            call split_path_name(file_base, path, file_base)
            l2auxFileHandle = GetPCFromRef(file_base, mlspcf_l2dgm_start, &
              & mlspcf_l2dgm_end, &
              & TOOLKIT, returnStatus, l2aux_Version, DEBUG, &
              & exactName=l2auxPhysicalFilename)
          else
            l2auxPhysicalFilename = file_base
            returnStatus = 0
          end if

          if ( returnStatus == 0 .and. .not. checkPaths ) then

            if ( DEBUG ) call output ( 'file name: ' // TRIM(l2auxPhysicalFilename), &
              & advance='yes' )
            ! Create the HDF file and initialize the SD interface
            if ( DEBUG ) call output ( 'Attempting sfstart', advance='yes' )
            sdfId = mls_sfstart(l2auxPhysicalFilename, DFACC_CREATE, &
              & hdfVersion=hdfVersion)

            if ( DEBUG ) call output ( "looping over quantities", advance='yes' )
            numquantitiesperfile = 0
            do field_no = 2, nsons(key) ! Skip "output" name
              gson = subtree(field_no,key)
              select case ( decoration(subtree(1,gson)) )
              case ( f_quantities )
                do in_field_no = 2, nsons(gson)
                  if ( DEBUG ) &
                    & call output ( "computing db index", advance='yes')
                  db_index = -decoration(decoration(subtree(in_field_no ,gson)))
                  if ( db_index >= 1 ) then
                    call WriteL2AUXData ( l2auxDatabase(db_index), sdfid, returnStatus,&
                      & WriteCounterMAF = &
                      &   (writeCounterMAF .and. numquantitiesperfile == 0), &
                      & hdfVersion=hdfVersion )
                    error = max(error, returnStatus)
                    numquantitiesperfile = numquantitiesperfile+1
                    if ( DEBUG ) call output(&
                      & "attempting to fill quantity name", advance='yes')
                    if ( numquantitiesperfile > MAXQUANTITIESPERFILE ) then
                      call announce_error ( son, &
                        & 'Attempt to write too many l2aux quantities to a file', &
                        & numquantitiesperfile )
                      numquantitiesperfile = MAXQUANTITIESPERFILE
                    end if
                    call get_string &
                      & ( l2auxDatabase(db_index)%name, &
                      &     QuantityNames(numquantitiesperfile) )
                  else
                    call MLSMessage ( MLSMSG_Warning, ModuleName, &
                      & 'Unable to save l2aux quantity, perhaps no chunks processed' )
                  end if
                end do ! in_field_no = 2, nsons(gson)
              case ( f_overlaps )
                ! ??? More work needed here
              end select
            end do ! field_no = 2, nsons(key)

            ! Now close the file
            returnStatus = mls_sfend(sdfid, hdfVersion=hdfVersion)
            !        returnStatus = sfend(sdfid)

            if ( returnStatus /= PGS_S_SUCCESS ) then
              call announce_error ( son, &
                &  "Error closing l2aux file:  "//l2auxPhysicalFilename, returnStatus)
            else if (index(switches, 'pro') /= 0) then
              call announce_success(l2auxPhysicalFilename, 'l2aux', &
                & numquantitiesperfile, quantityNames, hdfVersion=hdfVersion)
            end if

            if ( .not. TOOLKIT ) cycle

            ! Write the metadata file
            call add_metadata ( son, file_base, numquantitiesperfile, &
              & quantityNames, hdfVersion, l_hdf, metadata_error )

          else if ( returnStatus /= PGS_S_SUCCESS ) then
            call announce_error ( son, &
              &  "Error finding l2aux file matching:  "//file_base, returnStatus)
          end if

        case ( l_l2pc ) ! ------------------------------ Writing l2pc files --
          ! I intend to completely ignore the PCF file in this case,
          ! it's not worth the effort!
          ! In that case, I will ignore the possibility of checkPaths being true
          if ( .not. canWriteL2PC ) call MLSMessage(MLSMSG_Error,ModuleName,&
            & "Cannot write l2pc files with multi chunk l2cf's")
          recLen = 0
          packed = .false.
          ascii = .false.
          do field_no = 2, nsons(key) ! Skip "output" name
            gson = subtree(field_no,key)
            select case ( decoration(subtree(1,gson)) )
            case ( f_quantities )
              quantitiesNode = gson
            case ( f_overlaps )
              ! ??? More work needed here
            case ( f_packed )
              packed = get_boolean ( gson )
            case ( f_dontPack )
              call Allocate_Test ( dontPack, nsons(gson)-1, 'dontPack', ModuleName )
              do node = 2, nsons(gson)
                dontPack(node-1) = decoration(decoration(subtree(node,gson)))
              end do
            case ( f_ascii )
              ascii = get_boolean ( gson )
            end select
          end do ! field_no = 2, nsons(key)

          ! Open file
          if ( ascii ) then
            ! ASCII l2pc file
            l2pcUnit = mls_io_gen_openf ( 'open', .true., error,&
              & recLen, PGSd_IO_Gen_WSeqFrm, trim(file_base), 0,0,0, unknown=.true. )
            if ( error /= 0 ) call MLSMessage(MLSMSG_Error,ModuleName,&
              & 'Failed to open l2pc file:'//trim(file_base))

            do in_field_no = 2, nsons(quantitiesNode)
              db_index = decoration(decoration(subtree(in_field_no, quantitiesNode )))
              call GetFromMatrixDatabase ( matrices(db_index), tmpMatrix )
              call writeOneL2PC ( tmpMatrix, l2pcUnit, packed )
            end do ! in_field_no = 2, nsons(gson)

            error = mls_io_gen_closef ( 'cl', l2pcUnit)
            if ( error /= 0 ) then
              call MLSMessage(MLSMSG_Error,ModuleName,&
                & 'Failed to close l2pc file:'//trim(file_base))
            else if ( index(switches, 'pro') /= 0) then
              call announce_success(file_base, 'l2pc', &
                & 0, quantityNames)
            end if
          else
            ! For the moment call a routine
            call OutputHDF5L2PC ( trim(file_base), matrices, quantitiesNode, packed, &
              & dontPack )
            ! Later on when HDF5 is 'blessed' I want to move all this code
            ! here instead
            !            call H5FCreate_F ( trim(file_base), H5F_ACC_TRUNC, l2pcUnit, &
            !              & returnStatus )
            !            if ( returnStatus /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            !              & 'Unable to open hdf5 l2pc file for output.' )
            !            do in_field_no = 2, nsons(quantitiesNode)
            !              db_index = decoration(decoration(subtree(in_field_no, quantitiesNode )))
            !              call GetFromMatrixDatabase ( matrices(db_index), tmpMatrix )
            !              call writeOneHDF5L2PC ( tmpMatrix, l2pcUnit, packed )
            !            end do ! in_field_no = 2, nsons(gson)
            !            call H5FClose ( l2pcUnit, returnStatus )
            !            if ( returnStatus /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
            !              & 'Unable to close hdf5 l2pc file.' )

          end if
          call Deallocate_test ( dontPack, 'dontPack', ModuleName )

        case ( l_l2dgg ) ! --------------------- Writing l2dgg files -----

          if ( DEBUG ) call output('output file type l2dgg', advance='yes')
          ! Get the l2gp file name from the PCF

          if ( TOOLKIT ) then
            call split_path_name(file_base, path, file_base)
            l2gpFileHandle = GetPCFromRef(file_base, mlspcf_l2dgg_start, &
              & mlspcf_l2dgg_end, &
              & TOOLKIT, returnStatus, l2gp_Version, DEBUG, &
              & exactName=l2gpPhysicalFilename)
          else
            l2gpPhysicalFilename = file_base
            returnStatus = 0
          end if

         if ( returnStatus == 0 .and. .not. checkPaths ) then
            if ( DEBUG ) call output(&
              & 'file name: ' // TRIM(l2gpPhysicalFilename), advance='yes')
            ! Open the HDF-EOS file and write swath data

            if ( DEBUG ) call output('Attempting swopen', advance='yes')
            !           swfid = swopen(l2gpPhysicalFilename, DFACC_CREATE)
            swfid = mls_io_gen_openF('swopen', .TRUE., returnStatus, &
              & record_length, DFACC_CREATE, FileName=l2gpPhysicalFilename, &
              & hdfVersion=hdfVersion, debugOption=.false. )

            ! Loop over the segments of the l2cf line

            numquantitiesperfile = 0
            do field_no = 2, nsons(key) ! Skip "output" name
              gson = subtree(field_no,key)
              select case ( decoration(subtree(1,gson)) )
              case ( f_quantities )
                do in_field_no = 2, nsons(gson)
                  db_index = -decoration(decoration(subtree(in_field_no ,gson)))
                  call writeL2GPData ( l2gpDatabase(db_index), swfid, &
                    & hdfVersion=hdfVersion )
                  numquantitiesperfile = numquantitiesperfile+1
                  if ( numquantitiesperfile > MAXQUANTITIESPERFILE ) then
                    call announce_error ( son, &
                      & 'Attempt to write too many l2dgg quantities to a file', &
                      & numquantitiesperfile )
                    numquantitiesperfile = MAXQUANTITIESPERFILE
                  end if
                  quantityNames(numquantitiesperfile) = l2gpDatabase(db_index)%name
                end do ! in_field_no = 2, nsons(gson)
              case ( f_overlaps )
                ! ??? More work needed here
              end select
            end do ! field_no = 2, nsons(key)

            if ( DEBUG ) call output('Attempting swclose', advance='yes')
            !           returnStatus = swclose(swfid)
            returnStatus = mls_io_gen_closeF('swclose', swfid, &
              & hdfVersion=hdfVersion)
            if ( returnStatus /= PGS_S_SUCCESS ) then
              call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
              call MLSMessage ( MLSMSG_Error, ModuleName, &
                &  "Error closing  l2dgg file:  "//mnemonic//" "//msg )
            else if (index(switches, 'pro') /= 0) then
              call announce_success(l2gpPhysicalFilename, 'l2dgg', &
                & numquantitiesperfile, quantityNames, hdfVersion=hdfVersion)
            end if

            if ( .not. TOOLKIT ) cycle

            ! Write the metadata file
            call add_metadata ( son, file_base, numquantitiesperfile, &
              & quantityNames, hdfVersion, output_type, metadata_error )

          else if ( returnStatus /= PGS_S_SUCCESS ) then
            call announce_error ( son, &
              &  "Error finding l2gp file matching:  "//file_base, returnStatus)
          end if


        case default
          if ( any(output_type /= (/ l_l2gp, l_l2dgg, l_l2aux, l_l2pc /)) ) then
            call announce_error ( son, &
            &  "Error--unknown output type: parser should have caught this")
          else
            call output('Lahey did weird thing again: ', advance='yes')
          end if
            call output('l2gp type number: ', advance='no')
            call output(l_l2gp, advance='yes')

            call output('l2aux type number: ', advance='no')
            call output(l_l2aux, advance='yes')

            call output('l2dgg type number: ', advance='no')
            call output(l_l2dgg, advance='yes')

            call output('output type number: ', advance='no')
            call output(output_type, advance='yes')

            call output('file_base: ', advance='no')
            call output(trim(file_base), advance='yes')

        end select

      case ( s_time )
        if ( timing ) then
          call sayTime
        else
          call time_now ( t1 )
          timing = .true.
        end if

      case default
        call announce_error ( spec_no, &
          &  "Error--unknown spec_no: parser should have caught this")

      end select

    end do  ! spec_no
    
    ! Catenate any split Direct Writes
    ! We assume hdfVersion is 5
    if ( CATENATESPLITS .and. associated(DirectDatabase) ) then
    !! if ( .true. .and. associated(DirectDatabase) ) then
      ! Any dgg eligible for being catenated
      DB_index = findFirst( DirectDatabase%autoType, l_l2dgg )
      if ( findNext(DirectDatabase%autoType, l_l2dgg, DB_index) > 0 ) then
        if ( TOOLKIT ) then
          l2gp_Version = 1
          l2gpFileHandle = GetPCFromRef('DGG', mlspcf_l2dgg_start, &
            & mlspcf_l2dgg_end, &
            & TOOLKIT, returnStatus, l2gp_Version, DEBUG, &
            & exactName=l2gpPhysicalFilename)
        else
          file_base = DirectDatabase(DB_index)%fileName
          l2gpPhysicalFilename = unSplitName(file_base)
          returnStatus = 0
        end if
        if ( any(DirectDatabase%fileName == l2gpPhysicalFilename) ) then
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Cannot unsplit dgg dw to existing file " // &
            & trim(l2gpPhysicalFilename) )
        end if
        madeFile = .false.
        do DB_index = 1, size(DirectDatabase)
          if ( DirectDatabase(DB_index)%autoType /= l_l2dgg ) cycle
          if ( DEBUG ) then
            call output ( 'preparing to cp split dgg', advance='yes' )
            call output ( 'from: ', advance='no' )
            call output ( trim(DirectDatabase(DB_index)%fileName) , advance='yes' )
            call output ( '   to: ', advance='no' )
            call output ( trim(l2gpPhysicalFilename) , advance='yes' )
          end if
          if ( mls_exists(trim(DirectDatabase(DB_index)%fileName)) /= 0 ) cycle
          madeFile = .true.
          call cpL2GPData(trim(DirectDatabase(DB_index)%fileName), &
            & trim(l2gpPhysicalFilename), create2=(DB_index==1), &
            & hdfVersion1=HDFVERSION_5, hdfVersion2=HDFVERSION_5, &
            & notUnlimited=avoidUnlimitedDims, ReadStatus=.true.)
        end do
        if ( TOOLKIT .and. madeFile ) then
          call add_metadata ( 0, trim(l2gpPhysicalFilename), 1, &
            & (/'dgg'/), HDFVERSION_5, l_l2dgg, returnStatus )
          if ( returnStatus /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'unable to addmetadata to ' // trim(l2gpPhysicalFilename) )
        end if
      end if
      ! Next we would do the same for any split dgm direct write files
      DB_index = findFirst( DirectDatabase%autoType, l_l2aux )
      if ( findNext(DirectDatabase%autoType, l_l2aux, DB_index) > 0 ) then
        if ( TOOLKIT ) then
          l2gp_Version = 1
          l2gpFileHandle = GetPCFromRef('DGM', mlspcf_l2dgm_start, &
            & mlspcf_l2dgm_end, &
            & TOOLKIT, returnStatus, l2gp_Version, DEBUG, &
            & exactName=l2auxPhysicalFilename)
        else
          file_base = DirectDatabase(DB_index)%fileName
          l2auxPhysicalFilename = unSplitName(file_base)
          returnStatus = 0
        end if
        if ( any(DirectDatabase%fileName == l2auxPhysicalFilename) ) then
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            &  "Must not unsplit dgm dw to " // trim(l2auxPhysicalFilename) )
        end if
        madeFile = .false.
        do DB_index = 1, size(DirectDatabase)
          if ( DirectDatabase(DB_index)%autoType /= l_l2aux ) cycle
          ! print *, 'About to try to convert array2List'
          ! call dump(DirectDatabase(DB_index))
          call Array2List(DirectDatabase(DB_index)%sdNames, sdList)
          ! print *, 'result: ', trim(sdList)
          ! Not implemented yet
          if ( DEBUG ) then
            call output ( 'preparing to cp split dgm', advance='yes' )
            call output ( 'from: ', advance='no' )
            call output ( trim(DirectDatabase(DB_index)%fileName) , advance='yes' )
            call output ( '   to: ', advance='no' )
            call output ( trim(l2auxPhysicalFilename) , advance='yes' )
            call output ( '   sdList: ', advance='no' )
            call output ( trim(sdList) , advance='yes' )
          end if
          if ( mls_exists(trim(DirectDatabase(DB_index)%fileName)) /= 0 ) cycle
          madeFile = .true.
          if ( sdList /= ' ' ) then
            call cpL2AUXData(trim(DirectDatabase(DB_index)%fileName), &
            & trim(l2auxPhysicalFilename), create2=(DB_index==1), &
            & hdfVersion=HDFVERSION_5, sdList=trim(sdList))
          else
            ! Last-ditch effort if somehow sdNames empty or Array2List fails
            call cpL2AUXData(trim(DirectDatabase(DB_index)%fileName), &
            & trim(l2auxPhysicalFilename), create2=(DB_index==1), &
            & hdfVersion=HDFVERSION_5)
          end if
        end do
        ! Is metadata really needed for l2aux files?
        if ( TOOLKIT .and. madeFile ) then
          call add_metadata ( 0, trim(l2auxPhysicalFilename), 1, &
            & (/'dgm'/), HDFVERSION_5, l_hdf, returnStatus )
          if ( returnStatus /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'unable to addmetadata to ' // trim(l2auxPhysicalFilename) )
        end if
      end if
    end if

    ! Write the log file metadata
    if ( LOGFILEGETSMETADATA .and. .not. checkPaths ) then
      if ( DEBUG ) then
        call output('About to write log file metadata' , advance='yes')
      end if

      if ( TOOLKIT ) then
        call writeMetaLog ( metadata_error )
        error = max(error, PENALTY_FOR_NO_METADATA*metadata_error)
      end if
    end if

    ! Done with text of PCF file at last
   
    if ( DEBUG ) &
      & call output ( 'About to deallocate text of PCF file' , advance='yes' )

    call deallocate_test ( l2pcf%anText, 'anText of PCF file', moduleName )

    if (index(switches, 'pro') /= 0) then
      call output ( '============ End Level 2 Products ============', advance='yes' )
      call output ( ' ', advance='yes' )
    end if

    if ( error /= 0 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Problem with Output_Close section' )
    end if

    if ( timing ) call sayTime
    if ( toggle(gen) ) call trace_end ( "Output_Close")

  contains
    subroutine SayTime
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      end if
      call output ( "Timing for Output_Close =" )
      call output ( DBLE(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

  end subroutine Output_Close

  ! ---------------------------------------------  add_metadata  -----
  subroutine add_metadata ( node, fileName, numquantitiesperfile, &
    & quantityNames, hdfVersion, filetype, metadata_error )
    
    use INIT_TABLES_MODULE, only: L_L2DGG
    use Intrinsic, only: l_swath, l_hdf
    use MLSFiles, only: GetPCFromRef, split_path_name
    use MLSPCF2, only: MLSPCF_L2DGM_END, MLSPCF_L2DGM_START, MLSPCF_L2GP_END, &
      & MLSPCF_L2GP_START, mlspcf_l2dgg_start, mlspcf_l2dgg_end, &
      & Mlspcf_mcf_l2dgm_start, Mlspcf_mcf_l2dgg_start
    use WriteMetadata, only: Populate_metadata_std, &
      & Populate_metadata_oth, Get_l2gp_mcf
  ! Deal with metadata--1st for direct write, but later for all cases
  integer, intent(in) :: node
  character(len=*), intent(in) :: fileName
  integer, intent(in) :: numquantitiesperfile
  character(len=*), dimension(:), intent(in) :: quantityNames
  integer, intent(in) :: hdfVersion
  integer, intent(in) :: fileType  ! l_swath, l_hdf, ..
  integer, intent(out) :: metadata_error
  
  ! Internal variables
  integer :: baseIndex
  logical, parameter :: DEBUG = .false.
  character(len=*), parameter :: L2GPHEAD = 'L2GP-'
  integer :: field_no
  character (len=132) :: FILE_BASE
  integer :: fileHandle
  integer :: L2aux_mcf
  integer :: l2dgg_mcf
  integer :: L2gp_mcf
  character (len=32) :: meta_name=' '
  character (len=132) :: path
  character(len=len(fileName)) :: PhysicalFilename
  integer :: returnStatus
  integer :: Version
  ! Executable
  l2aux_mcf = mlspcf_mcf_l2dgm_start
  l2dgg_mcf = mlspcf_mcf_l2dgg_start
  metadata_error = 0
  Version = 1
  select case (filetype)
  case (l_swath, l_l2dgg)
     call split_path_name(fileName, path, file_base)
     baseIndex = index(trim(file_base), L2GPHEAD)
     if ( baseIndex > 0 ) then
       file_base=file_base(baseIndex+len(L2GPHEAD):)
     end if
     if ( filetype == l_l2dgg ) then
       FileHandle = GetPCFromRef(file_base, mlspcf_l2dgg_start, &
         & mlspcf_l2dgg_end, &
         & .true., returnStatus, Version, DEBUG, &
         & exactName=PhysicalFilename)
       l2gp_mcf = l2dgg_mcf
     else
       FileHandle = GetPCFromRef(file_base, mlspcf_l2gp_start, &
         & mlspcf_l2gp_end, &
         & .true., returnStatus, Version, DEBUG, &
         & exactName=PhysicalFilename)
       call get_l2gp_mcf ( file_base, meta_name, l2gp_mcf  )
     end if
     if (returnStatus /= 0) then
         call MLSMessage ( MLSMSG_Error, ModuleName, &
           &  "While adding metadata failed to GetPCFromRef for " // trim(fileName) )
     else if ( l2gp_mcf <= -999 ) then
         call MLSMessage ( MLSMSG_Warning, ModuleName, &
           &  "no mcf for this l2gp species in" // trim(file_base) )
         return
     else if (l2gp_mcf <= 0) then
         call MLSMessage ( MLSMSG_Error, ModuleName, &
           &  "no mcf for this l2gp species in" // trim(file_base) )
     end if

     if ( QuantityNames(numquantitiesperfile) &
       & == QuantityNames(1) ) then
       ! Typical homogeneous l2gp file: 
       ! e.g., associated with BrO is ML2BRO.001.MCF
       if ( DEBUG ) then
         call output('preparing to populate metadata_std', advance='yes')
         call output('l2gpFileHandle: ', advance='no')
         call output(FileHandle , advance='no')
         call output('   l2gp_mcf: ', advance='no')
         call output(l2gp_mcf , advance='yes')
       end if
       call populate_metadata_std &
         & (FileHandle, l2gp_mcf, QuantityNames(1), &
         & hdfVersion=hdfVersion, metadata_error=metadata_error, &
         & filetype=filetype  )
     else
       ! Type l2gp file 'other'
       if ( DEBUG ) then
         call output ( 'preparing to populate metadata_oth', advance='yes' )
         call output ( 'l2gpFileHandle: ', advance='no' )
         call output ( FileHandle , advance='no' )
         call output ( '   l2gp_mcf: ', advance='no' )
         call output ( l2gp_mcf , advance='yes' )
       end if

       call populate_metadata_oth &
         & ( FileHandle, l2gp_mcf, &
         & numquantitiesperfile, QuantityNames, &
         & hdfVersion=hdfVersion, metadata_error=metadata_error, &
         & filetype=filetype  )
     end if
  case (l_hdf)
     if ( DEBUG ) call output ( 'output file type l2aux', advance='yes' )
     ! Get the l2aux file name from the PCF

     call split_path_name(fileName, path, file_base)
     FileHandle = GetPCFromRef(file_base, mlspcf_l2dgm_start, &
       & mlspcf_l2dgm_end, &
       & .true., returnStatus, Version, DEBUG, &
       & exactName=PhysicalFilename)
     if (returnStatus /= 0) then
         call MLSMessage ( MLSMSG_Error, ModuleName, &
           &  "While adding metadata failed to GetPCFromRef for " // trim(fileName) )
     else if ( DEBUG ) then
       call output ( 'preparing to populate metadata_oth', advance='yes' )
       call output ( 'l2auxFileHandle: ', advance='no' )
       call output ( FileHandle , advance='no' )
       call output ( '   l2aux_mcf: ', advance='no' )
       call output ( l2aux_mcf , advance='no' )
       call output ( '   number of quantities: ', advance='no' )
       call output ( numquantitiesperfile , advance='yes' )
       do field_no=1, numquantitiesperfile
         call output ( field_no , advance='no' )
         call output ( '       ', advance='no' )
         call output ( trim(QuantityNames(field_no)) , advance='yes' )
       end do
     end if
     call populate_metadata_oth &
       & ( FileHandle, l2aux_mcf, &
       & numquantitiesperfile, QuantityNames,&
       & hdfVersion=hdfVersion, metadata_error=metadata_error, &
       & filetype=filetype  )
  case default
    call announce_error ( node, &
      &  "Error--filetype unrecognized", filetype)
    call MLSMessage(MLSMSG_Error,ModuleName,&
      & "Unrecognized filetype in add_metadata (must be swath or hdf) "&
      & // trim(filename))
  end select

  end subroutine add_metadata

! =====     Private Procedures     =====================================

  ! ---------------------------------------------  announce_success  -----
  subroutine announce_success ( Name, l2_type, num_quants, quantities, hdfVersion )
    integer, intent(in) :: num_quants 
    character(LEN=*), intent(in)   :: Name
    character(LEN=*), intent(in)   :: l2_type
    integer, optional,  intent(in) :: hdfVersion
    character(LEN=*), dimension(:), intent(in) :: quantities
    integer :: i

    call output ( 'Level 2 output product type : ' )
    call output ( trim(l2_type), advance='no')
    if ( present(hdfVersion) ) then
      call blanks(4)
      call output ( 'hdf ' )
      call output ( hdfVersion, advance='yes')
    else
      call output ( ' ', advance='yes')
    end if
    call blanks(15)
    call output ( 'name : ' )
    call blanks(8)
    call output ( trim(Name), advance='yes')

    if ( num_quants > 0 ) then
      call output ( 'number ' )
      call blanks(5)
      call output ( 'quantity', advance='yes')
      do i=1, num_quants
          call output ( i )
          call blanks(5)
          call output ( trim(quantities(i)), advance='yes')
      end do
    end  if
  end subroutine announce_success

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( Where, Full_message, Code, Penalty )

    use LEXER_CORE, only: PRINT_SOURCE
    use TREE, only: SOURCE_REF

    integer, intent(in) :: Where   ! Tree node where error was noticed
    character(LEN=*), intent(in) :: Full_Message
    integer, intent(in), optional :: Code    ! Code for error message
    integer, intent(in), optional :: Penalty
    integer :: myPenalty

    myPenalty = 1
    if ( present(penalty) ) myPenalty = penalty
    error = max(error,myPenalty)

    call output ( '***** At ' )
    if ( where > 0 ) then
      call print_source ( source_ref(where) )
    else
      call output ( '(no lcf node available)' )
    end if
    call output ( ' OutputAndClose complained: ' )

    call output ( trim(full_message), advance='yes', &
      & from_where=ModuleName )
    if ( present(code) ) then
      ! select case ( code )
      ! end select
      call output ( ' Code: ' )
      call output ( Code, advance='yes' )
    end  if
  end subroutine ANNOUNCE_ERROR

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module OutputAndClose

! $Log$
! Revision 2.98  2004/06/10 00:58:45  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.97  2004/05/19 20:22:09  vsnyder
! Remove USEs for unreferenced symbols, polish some cannonballs
!
! Revision 2.96  2004/04/24 00:22:55  pwagner
! Forcibly sets ReadStatus when knitting DGG; prunes L2GP- prefix from file_base
!
! Revision 2.95  2004/03/12 00:39:37  pwagner
! Interface to cpL2GPData changed to match
!
! Revision 2.94  2004/03/03 19:23:48  pwagner
! Master task never knows actual sdList to catenate; let cpL2AUXData figure it out
!
! Revision 2.93  2004/02/19 23:57:47  pwagner
! Hopefully will not try to write metadata if no file exists
!
! Revision 2.92  2004/02/05 23:35:21  pwagner
! Fixed some bugs in catenating split dgg/dgm directwrites
!
! Revision 2.91  2004/01/27 21:37:26  pwagner
! Can catenate split l2aux files
!
! Revision 2.90  2004/01/23 01:15:00  pwagner
! Began effort to catenate split direct write files
!
! Revision 2.89  2004/01/22 00:52:19  pwagner
! Small changes regarding metadata
!
! Revision 2.88  2003/11/07 00:46:51  pwagner
! New quicker preflight option: --checkPaths
!
! Revision 2.87  2003/10/28 21:42:36  pwagner
! Exits with message if cant GetPCFFromRef
!
! Revision 2.86  2003/10/20 23:59:20  pwagner
! Simplified code for writing metadata
!
! Revision 2.85  2003/10/16 18:29:35  pwagner
! Should not try to write metadata twice on DirectWrite files
!
! Revision 2.84  2003/09/19 23:29:27  pwagner
! Should not be a metadata error when DirectWrite-ing CH3CN
!
! Revision 2.83  2003/09/04 22:42:47  pwagner
! Some tweaks relating to DirectWrite; may not matter
!
! Revision 2.82  2003/08/15 20:43:10  pwagner
! Downgraded to warning if directwrite output_type unkown, e.g. l_l2fwm
!
! Revision 2.81  2003/08/08 23:06:39  livesey
! Added the dontPack option on outputing l2pc files.
!
! Revision 2.80  2003/08/01 20:07:44  pwagner
! Fixed Toolkit bug; metadata distinguishes l2dgg from l2gp
!
! Revision 2.79  2003/07/23 18:29:32  cvuu
! quick and dirty fixed for CH3CN
!
! Revision 2.78  2003/07/07 23:49:11  pwagner
! Add_metadata procedure now public
!
! Revision 2.77  2003/07/07 17:31:44  livesey
! Mainly cosmetic changes
!
! Revision 2.76  2003/06/26 00:17:17  pwagner
! Writes metadata to all files in DirectDataBase
!
! Revision 2.75  2003/06/24 23:54:07  pwagner
! New db indexes stored for entire direct file
!
! Revision 2.74  2003/06/23 18:06:33  pwagner
! Should allow us to write metadata after DirectWrite
!
! Revision 2.73  2003/06/20 19:38:26  pwagner
! Allows direct writing of output products
!
! Revision 2.72  2003/06/09 22:49:33  pwagner
! Reduced everything (PCF, PUNISH.., etc.) to TOOLKIT
!
! Revision 2.71  2003/05/12 02:07:06  livesey
! Bound r8->r4 conversion in direct write
!
! Revision 2.70  2003/04/03 22:59:23  pwagner
! setAlias no longer an arg to write_meta
!
! Revision 2.69  2003/03/20 19:22:56  pwagner
! Fixed bug in DirectWrite_hdf5; seems to work
!
! Revision 2.68  2003/03/11 00:21:36  pwagner
! Interfaces fit new WritePCF2Hdr flixibility
!
! Revision 2.67  2003/03/01 00:25:20  pwagner
! Disabled writing metadata to Log file (aka PH)
!
! Revision 2.66  2003/02/12 21:51:32  pwagner
! Should allow direct write with attributes
!
! Revision 2.65  2003/02/10 22:01:54  pwagner
! Passes isHDFEOS to metadata; writes globalattributes during DirectWrite
!
! Revision 2.64  2003/01/23 23:31:42  pwagner
! May directwrite to hdf5 l2aux files
!
! Revision 2.63  2002/12/11 22:21:05  pwagner
! Makes soft link to data field name from L2gpValue field in hdf5 l2gp
!
! Revision 2.62  2002/11/22 19:10:30  pwagner
! Upped MAXQUANTITIESPERFILE to 10k
!
! Revision 2.61  2002/11/13 01:10:09  pwagner
! Beginnings of attempt to write hdf5 L2AUX; incomplete
!
! Revision 2.60  2002/10/08 17:36:22  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.59  2002/08/21 02:35:18  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.58  2002/08/21 01:05:06  livesey
! Changed to single precision for direct write
!
! Revision 2.57  2002/08/20 04:37:06  livesey
! Minor typo
!
! Revision 2.56  2002/08/20 04:33:13  livesey
! Added extra check in direct write
!
! Revision 2.55  2002/08/15 21:47:04  pwagner
! WriteL2AuxData now returns non-zero status if it fails
!
! Revision 2.54  2002/06/12 17:58:42  livesey
! Intermediate support for HDF5 L2PCs
!
! Revision 2.53  2002/05/22 16:30:31  livesey
! Bug fix in directWrite
!
! Revision 2.52  2002/05/22 00:49:01  livesey
! Added direct write stuff
!
! Revision 2.51  2002/05/07 20:26:15  livesey
! Added writeCounterMAF option for l2aux
!
! Revision 2.50  2002/02/22 19:19:48  pwagner
! Fixed bug in metaName use
!
! Revision 2.49  2002/02/22 01:16:17  pwagner
! Uses new metaName field for mcf file hint
!
! Revision 2.48  2002/01/29 23:49:38  pwagner
! Separate DEFAULT_HDFVERSION_(READ)(WRITE)
!
! Revision 2.47  2002/01/26 00:10:45  pwagner
! Correctly sets hdfVersion; changed proclaim to announce_success
!
! Revision 2.46  2002/01/23 21:52:15  pwagner
! Accepts and uses hdfVersion optional field
!
! Revision 2.45  2002/01/18 23:07:48  pwagner
! Uses MLSFiles instead of HDFEOS
!
! Revision 2.44  2002/01/18 00:24:34  livesey
! Added packed option to outputing l2pc files
!
! Revision 2.43  2001/11/20 00:48:54  livesey
! Alleviated one bug in zero chunks case, but there's another one to
! fix later.  We need to decide how to handle this one.
!
! Revision 2.42  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.41  2001/11/01 21:05:10  pwagner
! Satisfies new WriteL2AuxData interface
!
! Revision 2.40  2001/10/12 23:12:23  pwagner
! Checks that number of quantities written to a file not grow too large
!
! Revision 2.39  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.38  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
! Revision 2.37  2001/06/04 23:57:40  pwagner
! Splits path from l2cf-defined file name before getPCfromRef
!
! Revision 2.36  2001/05/30 23:03:13  pwagner
! Moved PCFL2CFSAMECASE to MLSL2Options
!
! Revision 2.35  2001/05/17 22:33:28  pwagner
! Prints info if pro switch set
!
! Revision 2.34  2001/05/04 23:22:13  pwagner
! Detachable from Toolkit; created metafiles conditionally
!
! Revision 2.33  2001/05/04 23:19:55  pwagner
! Detachable from Toolkit; created metafiles conditionally
!
! Revision 2.32  2001/05/03 20:32:33  vsnyder
! Add a nullify and some cosmetic changes
!
! Revision 2.31  2001/05/01 23:57:23  pwagner
! Added l2dgg output type
!
! Revision 2.30  2001/04/28 01:30:52  livesey
! Stuff formerly outputting L2PCs is now outputting matrices.
!
! Revision 2.29  2001/04/26 20:02:09  livesey
! Made l2pc database a saved array in L2PC_m
!
! Revision 2.28  2001/04/26 15:59:13  livesey
! Fixed arguments to writeOneL2PC
!
! Revision 2.27  2001/04/25 21:51:28  livesey
! Minor changes, add canWriteL2PC flag
!
! Revision 2.26  2001/04/25 20:34:04  livesey
! Now writes l2pc files
!
! Revision 2.25  2001/04/20 23:51:24  vsnyder
! Improve an error message.  Add an option to consider the penalty in
! Announce_Error.  Numerous cosmetic changes.
!
! Revision 2.24  2001/04/19 23:51:40  pwagner
! Moved anText to become component of PCFData_T
!
! Revision 2.23  2001/04/16 23:51:08  pwagner
! Gets penalty from MLSL2Options
!
! Revision 2.22  2001/04/13 23:48:37  pwagner
! Removed bogus increment of l2gp_mcf
!
! Revision 2.21  2001/04/13 00:26:23  pwagner
! Whether files named in PCF agree in case with l2cf controlled
!
! Revision 2.19  2001/04/10 23:01:57  pwagner
! Now works better; tacks if no metadata
!
! Revision 2.18  2001/04/09 23:44:34  pwagner
! Fewer mistakes, more debug-type output
!
! Revision 2.17  2001/04/07 00:13:44  pwagner
! Extra error checks
!
! Revision 2.16  2001/04/04 23:44:52  pwagner
! Now deallocates anText of PCF file at last
!
! Revision 2.15  2001/04/03 23:51:28  pwagner
! Many changes; some may be right
!
! Revision 2.14  2001/04/02 23:43:46  pwagner
! Now makes metadata calls; it compiles, but does it bomb?
!
! Revision 2.13  2001/03/28 00:23:20  pwagner
! Made tiny changes to use announce_error
!
! Revision 2.12  2001/03/20 18:35:02  pwagner
! Using GetPCFromRef to get file handles
!
! Revision 2.11  2001/03/15 21:18:57  vsnyder
! Use Get_Spec_ID instead of decoration(subtree...
!
! Revision 2.10  2001/03/06 22:40:24  livesey
! Working l2aux
!
! Revision 2.9  2001/02/23 18:15:48  livesey
! Added trace calls.
!
! Revision 2.8  2001/02/13 22:59:36  pwagner
! l2 modules can only use MLSPCF2
!
! Revision 2.7  2001/02/09 19:30:16  vsnyder
! Move checking for required and duplicate fields to init_tables_module
!
! Revision 2.6  2001/02/09 00:38:22  livesey
! Various updates
!
! Revision 2.5  2001/01/03 18:15:13  pwagner
! Changed types of t1, t2 to real
!
! Revision 2.4  2000/11/16 02:25:13  vsnyder
! Implement timing.
!
! Revision 2.3  2000/10/05 16:43:00  pwagner
! Now compiles with new L2GPData module
!
! Revision 2.2  2000/09/11 19:43:47  ahanzel
! Removed old log entries.
!
! Revision 2.1  2000/09/08 22:55:56  vsnyder
! Revised to use the tree output by the parser
!

