! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=======================================================================================

module OutputAndClose ! outputs all data from the Join module to the
                      ! appropriate L2 Files

!=======================================================================================

  use Allocate_Deallocate, only: Deallocate_Test
  use Hdf, only: DFACC_CREATE, SFEND, SFSTART
  use HDFEOS, only: SWCLOSE, SWOPEN
  use INIT_TABLES_MODULE, only: F_FILE, F_OVERLAPS, F_QUANTITIES, F_TYPE, &
    & L_L2AUX, L_L2GP, L_L2PC, LIT_INDICES, S_OUTPUT, S_TIME
  use L2AUXData, only: L2AUXDATA_T, WriteL2AUXData
  use L2GPData, only: L2GPData_T, WriteL2GPData, L2GPNameLen
  use L2PC_m, only: L2PC_T, WRITEONEL2PC
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: I4
  use MLSFiles, only: GetPCFromRef, MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF
  use MLSL2Options, only: PENALTY_FOR_NO_METADATA
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSPCF2, only: MLSPCF_L2DGM_END, MLSPCF_L2DGM_START, MLSPCF_L2GP_END, &
    & MLSPCF_L2GP_START, &
    & Mlspcf_mcf_l2gp_start, Mlspcf_mcf_l2dgm_start, &
    & Mlspcf_mcf_l2dgg_start
  use MoreTree, only: Get_Spec_ID
  use OUTPUT_M, only: OUTPUT
  use SDPToolkit, only: PGS_S_SUCCESS, PGS_SMF_GETMSG, PGSD_IO_GEN_WSEQFRM
  use STRING_TABLE, only: GET_STRING
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TOGGLES, only: GEN, TOGGLE
  use TREE, only: DECORATION, DUMP_TREE_NODE, NODE_ID, NSONS, SOURCE_REF, &
    & SUBTREE, SUB_ROSA
  use TREE_TYPES, only: N_NAMED
  use WriteMetadata, only: PCFData_T, Populate_metadata_std, &
    & Populate_metadata_oth, WriteMetaLog, Get_l2gp_mcf

  implicit none
  private
  public :: Output_Close

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = &
    & "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  !-----------------------------------------------------------------------------

  ! -----     Private declarations     ---------------------------------

  ! For Announce_Error
  integer :: ERROR

  ! Must files named in PCF have same case as short names used in l2cf?
  logical, parameter :: PCFL2FCSAMECASE = .FALSE.

contains ! =====     Public Procedures     =============================

  ! -----------------------------------------------  Output_Close  -----
  subroutine Output_Close ( root, l2gpDatabase, l2auxDatabase, l2pcDatabase,&
    & l2pcf )

    ! Hard-wired assumptions:

    ! ----------------------- metadata ------------------------

    !   for the l2aux the mcf is mlspcf_mcf_l2dgm_start
    !   for the log file the mcf is mlspcf_mcf_l2log_start
    !   for the dgg file the mcf is mlspcf_mcf_l2dgg_start

    ! The correspondence between MCF and l2gp files is determined by
    ! the value of        MCFFORL2GPOPTION
    ! (see write_metadata module for fuller explanation)

    ! Arguments
    integer, intent(in) :: ROOT   ! Of the output section's AST
    type (L2GPData_T), dimension(:), pointer :: l2gpDatabase ! L2GP products
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase ! L2AUX products
    type (L2PC_T), dimension(:), pointer :: l2pcDatabase ! L2PC products
    type(PCFData_T) :: l2pcf

    ! - - - Local declarations - - -

    integer :: DB_index
    integer :: FIELD_INDEX              ! F_... field code
    integer :: FIELD_NO                 ! Index of assign vertex sons of Key
    character (len=132) :: FILE_BASE    ! From the FILE= field
    integer :: GSON                     ! Son of Son -- an assign node
    integer :: IN_FIELD_NO              ! Index of sons of assign vertex
    integer :: KEY                      ! Index of spec_args node
    integer :: l2auxFileHandle, l2aux_Version
    integer :: l2gp_mcf, l2aux_mcf, l2dgg_mcf  ! mcf numbers for writing metadata
    integer :: L2PCUNIT
    character (len=132) :: l2auxPhysicalFilename
    integer :: l2gpFileHandle, l2gp_Version
    character (len=132) :: l2gpPhysicalFilename
    integer, parameter:: MAXQUANTITIESPERFILE=64        
    integer :: metadata_error
    character (len=32) :: mnemonic
    character (len=256) :: msg
    integer :: NAME                     ! string index of label on output
    integer:: numquantitiesperfile        
    integer :: OUTPUT_TYPE              ! L_L2AUX or L_L2GP
    character(len=L2GPNameLen), dimension(MAXQUANTITIESPERFILE) :: QuantityNames  ! From "quantities" field
    integer :: RECLEN                   ! For file stuff
    integer :: returnStatus
    integer(i4) :: SDFID                ! File handle
    integer :: SON                      ! Of Root -- spec_args or named node
    integer :: SPEC_NO                  ! Index of son of Root
    integer :: SWFID
    real :: T1, T2     ! for timing
    logical :: TIMING
    logical, parameter :: DEBUG = .FALSE.

    ! Executable code
    timing = .false.

    if ( toggle(gen) ) call trace_begin ( "Output_Close", root)

    error = 0

    ! l2gp_mcf will be incremented in get_l2gp_mcf (if MCFFORL2GPOPTION == 1)
    l2gp_mcf = mlspcf_mcf_l2gp_start - 1   

    l2aux_mcf = mlspcf_mcf_l2dgm_start
    l2dgg_mcf = mlspcf_mcf_l2dgg_start

    ! Loop over the lines in the l2cf

    do spec_no = 2, nsons(root)-1 ! Skip name at begin and end of section

      l2gp_Version = 1
      l2aux_Version = 1

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
          field_index = decoration(subtree(1,gson))
          select case ( field_index )   ! Field name
          case ( f_file )
            call get_string ( sub_rosa(subtree(2,gson)), file_base )
            file_base = file_base(2:LEN_TRIM(file_base)-1) ! Parser includes quotes
          case ( f_type )
            output_type = decoration(subtree(2,gson))
          case default                  ! Everything else processed later
          end select
        end do

        select case ( output_type )
        case ( l_l2gp ) ! --------------------- Writing l2gp files -----
          if ( DEBUG ) call output('output file type l2gp', advance='yes')
          ! Get the l2gp file name from the PCF

          l2gpFileHandle = GetPCFromRef(file_base, mlspcf_l2gp_start, &
            & mlspcf_l2gp_end, &
            & PCFL2FCSAMECASE, returnStatus, l2gp_Version, DEBUG, &
            & exactName=l2gpPhysicalFilename)

          if ( returnStatus == 0 ) then
            if ( DEBUG ) call output(&
              & 'file name: ' // TRIM(l2gpPhysicalFilename), advance='yes')
            ! Open the HDF-EOS file and write swath data

            if ( DEBUG ) call output('Attempting swopen', advance='yes')
            swfid = swopen(l2gpPhysicalFilename, DFACC_CREATE)

            ! Loop over the segments of the l2cf line

            numquantitiesperfile = 0
            do field_no = 2, nsons(key) ! Skip "output" name
              gson = subtree(field_no,key)
              select case ( decoration(subtree(1,gson)) )
              case ( f_quantities )
                do in_field_no = 2, nsons(gson)
                  db_index = -decoration(decoration(subtree(in_field_no ,gson)))
                  call writeL2GPData ( l2gpDatabase(db_index), swfid )
                  numquantitiesperfile = numquantitiesperfile+1
                  quantityNames(numquantitiesperfile) = l2gpDatabase(db_index)%name
                end do ! in_field_no = 2, nsons(gson)
              case ( f_overlaps )
                ! ??? More work needed here
              end select
            end do ! field_no = 2, nsons(key)

            if ( DEBUG ) call output('Attempting swclose', advance='yes')
            returnStatus = swclose(swfid)
            if ( returnStatus /= PGS_S_SUCCESS ) then
              call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
              call MLSMessage ( MLSMSG_Error, ModuleName, &
                &  "Error closing  l2gp file:  "//mnemonic//" "//msg )
            end if

            ! Write the metadata file

            call get_l2gp_mcf ( file_base, l2gp_mcf, l2pcf )

            if ( l2gp_mcf <= 0 ) then

              ! Error in finding mcf number
              call announce_error ( son, &
                & 'No mcf numbers correspond to this l2gp file', l2gp_mcf, &
                & PENALTY_FOR_NO_METADATA )

            else if ( numquantitiesperfile <= 0 ) then

              ! Error in number of quantities
              call announce_error ( son, &
                & 'No quantities written for this l2gp file')

            else if ( QuantityNames(numquantitiesperfile) &
              & == QuantityNames(1) ) then

              ! Typical homogeneous l2gp file: 
              ! e.g., associated with BrO is ML2BRO.001.MCF
              if ( DEBUG ) then
                call output('preparing to populate metadata_std', advance='yes')
                call output('l2gpFileHandle: ', advance='no')
                call output(l2gpFileHandle , advance='no')
                call output('   l2gp_mcf: ', advance='no')
                call output(l2gp_mcf , advance='no')
                call output('   swfid: ', advance='no')
                call output(swfid , advance='yes')
              end if

              call populate_metadata_std ( &
                & l2gpFileHandle, l2gp_mcf, l2pcf, QuantityNames(1), &
                & metadata_error )
              error = max(error, PENALTY_FOR_NO_METADATA*metadata_error)

            else

              ! Type l2gp file 'other'
              if ( DEBUG ) then
                call output ( 'preparing to populate metadata_oth', advance='yes' )
                call output ( 'l2gpFileHandle: ', advance='no' )
                call output ( l2gpFileHandle , advance='no' )
                call output ( '   l2gp_mcf: ', advance='no' )
                call output ( l2gp_mcf , advance='no' )
                call output ( '   swfid: ', advance='no' )
                call output ( swfid , advance='yes' )
              end if

              call populate_metadata_oth ( &
                & l2gpFileHandle, l2gp_mcf, l2pcf, &
                & numquantitiesperfile, QuantityNames, metadata_error )
              error = max(error, PENALTY_FOR_NO_METADATA*metadata_error)
            end if

          else
            call announce_error ( ROOT, &
              &  "Error finding l2gp file matching:  "//file_base, returnStatus)
          end if

        case ( l_l2aux ) ! ------------------------------ Writing l2aux files ---

          if ( DEBUG ) call output ( 'output file type l2aux', advance='yes' )
          ! Get the l2aux file name from the PCF

          l2auxFileHandle = GetPCFromRef(file_base, mlspcf_l2dgm_start, &
            & mlspcf_l2dgm_end, &
            & PCFL2FCSAMECASE, returnStatus, l2aux_Version, DEBUG, &
            & exactName=l2auxPhysicalFilename)

          if ( returnStatus == 0 ) then

            if ( DEBUG ) call output ( 'file name: ' // TRIM(l2auxPhysicalFilename), &
              & advance='yes' )
            ! Create the HDF file and initialize the SD interface
            if ( DEBUG ) call output ( 'Attempting sfstart', advance='yes' )
            sdfId = sfstart(l2auxPhysicalFilename, DFACC_CREATE)

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
                  if ( DEBUG ) then
                    call output ( "db index:  ", advance='no')
                    call output ( db_index, advance='yes')
                    call output ( "db size:  ", advance='no')
                    call output ( size(l2auxDatabase), advance='yes')
                  end if
                  if ( db_index <= size(l2auxDatabase) ) then
                    call WriteL2AUXData ( l2auxDatabase(db_index),sdfid )
                    numquantitiesperfile = numquantitiesperfile+1
                    if ( DEBUG ) call output(&
                      & "attempting to fill quantity name", advance='yes')
                    call get_string &
                      & ( l2auxDatabase(db_index)%name, &
                      &     QuantityNames(numquantitiesperfile) )
                  else
                    call announce_error ( root, &
              	      &  "l2aux database smaller than db_index  ")
                    call output ( "db size:  ", advance='no')
                    call output ( size(l2auxDatabase), advance='yes')
                    call output ( "db index:  ", advance='no')
                    call output ( db_index, advance='yes')
                  end if
                end do ! in_field_no = 2, nsons(gson)
              case ( f_overlaps )
                ! ??? More work needed here
              end select
            end do ! field_no = 2, nsons(key)

            ! Now close the file
            returnStatus = sfend(sdfid)
            if ( returnStatus /= PGS_S_SUCCESS ) then
              call announce_error ( root, &
                &  "Error closing l2aux file:  "//l2auxPhysicalFilename, returnStatus)
            end if

            ! Write the metadata file
            if ( numquantitiesperfile <= 0 ) then
	      call announce_error ( son, &
	        & 'No quantities written for this l2aux file')
            else

              ! We may need to think more about this; until then reuse
              ! populate_metadata_oth, but with l2aux_mcf
              if ( DEBUG ) then
                call output ( 'preparing to populate metadata_oth', advance='yes' )
                call output ( 'l2auxFileHandle: ', advance='no' )
                call output ( l2auxFileHandle , advance='no' )
                call output ( '   l2aux_mcf: ', advance='no' )
                call output ( l2aux_mcf , advance='no' )
                call output ( '   sdfId: ', advance='no' )
                call output ( sdfId , advance='yes' )
                call output ( '   number of quantities: ', advance='no' )
                call output ( numquantitiesperfile , advance='yes' )
                do field_no=1, numquantitiesperfile
                  call output ( field_no , advance='no' )
                  call output ( '       ', advance='no' )
                  call output ( trim(QuantityNames(field_no)) , advance='yes' )
                end do
              end if
              call populate_metadata_oth ( &
                & l2auxFileHandle, l2aux_mcf, l2pcf, &
                & numquantitiesperfile, QuantityNames, metadata_error )
              error = max(error, PENALTY_FOR_NO_METADATA*metadata_error)
            end if

          else
            call announce_error ( root, &
              &  "Error finding l2aux file matching:  "//file_base, returnStatus)
          end if

        case ( l_l2pc ) ! ------------------------------ Writing l2pc files --
          ! I intend to completely ignore the PCF file in this case,
          ! it's not worth the effort!
          recLen = 0
          l2pcUnit = mls_io_gen_openf ( 'open', .true., error,&
            & recLen, PGSd_IO_Gen_WSeqFrm, trim(file_base), 0,0,0 )
          if ( error /= 0 ) call MLSMessage(MLSMSG_Error,ModuleName,&
            & 'Failed to open l2pc file:'//trim(file_base))
          do field_no = 2, nsons(key) ! Skip "output" name
            gson = subtree(field_no,key)
            select case ( decoration(subtree(1,gson)) )
            case ( f_quantities )
              do in_field_no = 2, nsons(gson)
                db_index = decoration(decoration(subtree(in_field_no ,gson)))
                call writeOneL2PC ( l2pcDatabase(db_index), l2pcUnit, lit_indices )
              end do ! in_field_no = 2, nsons(gson)
            case ( f_overlaps )
              ! ??? More work needed here
            end select
          end do ! field_no = 2, nsons(key)
          error = mls_io_gen_closef ( 'cl', l2pcUnit)
          if ( error /= 0 ) call MLSMessage(MLSMSG_Error,ModuleName,&
            & 'Failed to close l2pc file:'//trim(file_base))
        end select

      case ( s_time )
        if ( timing ) then
          call sayTime
        else
          call cpu_time ( t1 )
          timing = .true.
        end if
      end select
    end do  ! spec_no
    if ( timing ) call sayTime

! Write the log file metadata

    if ( DEBUG ) then
      call output('About to write log file metadata' , advance='yes')
    end if

    call writeMetaLog ( l2pcf, metadata_error )
    error = max(error, PENALTY_FOR_NO_METADATA*metadata_error)

! Done with text of PCF file at last

    if ( DEBUG ) &
      & call output ( 'About to deallocate text of PCF file' , advance='yes' )

    call deallocate_test ( l2pcf%anText, 'anText of PCF file', moduleName )

    if ( error /= 0 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Problem with Output_Close section' )
    end if

    if ( toggle(gen) ) call trace_end ( "Output_Close")

  contains
    subroutine SayTime
      call cpu_time ( t2 )
      call output ( "Timing for Output_Close =" )
      call output ( DBLE(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

  end subroutine Output_Close

! =====     Private Procedures     =====================================

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( Where, Full_message, Code, Penalty )
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
      select case ( code )
      end select
    end  if
  end subroutine ANNOUNCE_ERROR

end module OutputAndClose

! $Log$
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

