! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ReadAPriori

  use GriddedData, only: GriddedData_T, v_is_pressure, &
    & AddGriddedDataToDatabase, Dump
  use Hdf, only: DFACC_READ, SFSTART
  use Hdfeos, only: SWOPEN, SWCLOSE
  use INIT_TABLES_MODULE, only: F_FIELD, F_FILE, F_ORIGIN, F_SDNAME, F_SWATH, &
    & FIELD_FIRST, FIELD_LAST, L_CLIMATOLOGY, L_DAO, L_NCEP, S_GRIDDED, &
    & S_L2AUX, S_L2GP
  use L2AUXData, only: L2AUXData_T, AddL2AUXToDatabase, &
    &                  ReadL2AUXData, Dump
  use L2GPData, only: L2GPData_T, AddL2GPToDatabase, ReadL2GPData, Dump
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: FileNameLen
  use MLSFiles, only: GetPCFromRef, MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF, &
    & split_path_name, mls_InqSwath
  use MLSL2Options, only: PCF
  use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSPCF2, only: mlspcf_l2clim_start, mlspcf_l2clim_end
  use MoreTree, only: Get_Spec_ID
  use ncep_dao, only: READ_CLIMATOLOGY, ReadGriddedData
  use OUTPUT_M, only: BLANKS, OUTPUT
  use SDPToolkit, only: Pgs_pc_getReference, PGS_S_SUCCESS
  use String_Table, only: GET_STRING
  use Time_M, only: Time_Now
  use TOGGLES, only: GEN, SWITCHES, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, &
    &             SUB_ROSA, SUBTREE, DUMP_TREE_NODE, SOURCE_REF
  use TREE_TYPES, only: N_NAMED

  implicit none
  private
  public ::  read_apriori
  private ::  announce_error
  integer, private :: ERROR

  ! -----     Private declarations     ---------------------------------

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! --------------------------------------------------  read_apriori  -----
  ! Read a priori data from data files, be they l2gp, l2aux, climatology,
  ! NCEP, DAO etc.

  subroutine Read_apriori ( Root, L2GPDatabase, L2auxDatabase, GriddedDatabase)

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the Read a priori section in the AST
    type (l2gpdata_t), dimension(:), pointer :: L2GPDatabase
    type (L2AUXData_T), dimension(:), pointer :: L2auxDatabase
    type (GriddedData_T), dimension(:), pointer :: GriddedDatabase 

    ! Local Variables
    integer :: COMMAPOS                 ! For parsing string
    integer :: Details             ! How much info about the files to dump
    integer :: FIELD               ! Son of KEY, must be n_assign
    integer :: FIELDINDEX          ! Literal
    integer :: FieldName        ! sub-rosa index of name in field='name'
    character(len=FileNameLen) :: FIELDNAMESTRING ! actual literal clim. field
    integer :: FileHandle          ! fileHandle of a priori data file
    integer :: FileName            ! Sub-rosa index of name in file='name'
    character(len=FileNameLen) :: FileNameString   ! actual literal file name
    integer :: FileType            ! either s_l2gp or s_l2aux
    logical, dimension(field_first:field_last) :: GOT
    type (griddedData_T) :: GriddedData
    integer :: GriddedOrigin            ! From tree
    integer :: GridIndex           ! In the griddeddata database
    logical :: GotAlready               ! Do we need to reread this file?
    integer :: I, J                ! Loop indices for section, spec
    integer :: KEY                 ! Index of n_spec_args in the AST
    integer :: LastClimPCF
    integer :: LISTSIZE                 ! Size of string from SWInqSwath
    type (L2AUXData_T) :: L2AUX
    type (L2GPData_T) :: L2GP
    integer :: L2Index             ! In the l2gp or l2aux database
    integer :: L2Name              ! Sub-rosa index of L2[aux/gp] label
    integer :: NOSWATHS                 ! In an input file
    integer :: pcf_indx            ! loop index of climatology pcf numbers
    integer :: record_length
    integer :: ReturnStatus
    integer :: SON              ! Of root, an n_spec_args or a n_named
    integer :: SdName        ! sub-rosa index of name in sdName='name'
    character(len=FileNameLen) :: SDNAMESTRING ! actual literal sdName
    integer :: Sd_id
    integer :: SwathName        ! sub-rosa index of name in swath='name'
    character(len=FileNameLen) :: SWATHNAMESTRING ! actual literal swath name
    real :: T1, T2                      ! for timing
    logical :: TIMING
    integer :: Version

    character(len=2048) :: ALLSWATHNAMES ! Buffer to get info back.

    if ( toggle (gen) ) call trace_begin ( "read_apriori", root )

    timing = section_times
    if ( timing ) call time_now ( t1 )
    error = 0

    ! Will we be dumping info? To what level of detail?
    if ( index(switches, 'apr3') /= 0 ) then
      Details = 1
    elseif ( index(switches, 'apr2') /= 0 ) then
      Details = 0
    elseif ( index(switches, 'apr1') /= 0 ) then
      Details = -1
    else
      Details = -2
    endif
    if( index(switches, 'apr') /= 0 ) &     
    & call output ( '============ Read APriori ============', advance='yes' )    
    version = 1
    allswathnames = ' '
    lastClimPCF = mlspcf_l2clim_start - 1

    do i = 2, nsons(root)-1 ! Skip the section name at begin and end
      got = .false.
      son = subtree(i,root)
      if ( node_id(son) == n_named ) then ! Is spec labeled?
        key = subtree(2,son)
        l2Name = sub_rosa(subtree(1,son))
      else
        key = son
        l2Name = 0
      end if

      ! Node_id(key) is now n_spec_args.

      FileType = get_spec_id(key)

      ! Now parse file and field names
      fileName = 0
      swathName = 0
      do j = 2, nsons(key)
        field = subtree(j,key)
        fieldIndex = decoration(subtree(1,field))
        got(fieldIndex) = .true.
        select case ( fieldIndex )
        case ( f_file )
          fileName = sub_rosa(subtree(2,field))
        case ( f_swath )
          swathName = sub_rosa(subtree(2,field))
        case ( f_sdname )
          sdname = sub_rosa(subtree(2,field))
        case ( f_field )
          fieldName = sub_rosa(subtree(2,field))
        case ( f_origin )
          griddedOrigin = decoration(subtree(2,subtree(j,key)))
        end select
      end do

      if ( got(f_file) ) then
        call get_string ( FileName, fileNameString, strip=.true. )
      else
        fileNameString = ''
      end if
        
      select case ( FileType )
      case ( s_l2gp )
        if ( .not. got(f_file) ) &
          & call announce_error ( son, &
            & 'Filename name must be specified in read a priori' )
        swathNameString=''
        if ( got(f_swath) ) &
          & call get_string ( swathName, swathNameString, strip=.true. )

        ! If we didn't get a name get the first swath name in the file
        if ( len_trim(swathNameString) == 0 ) then
!
! (((((( This will have to be changed before transition to hdfeos5 ))))))
!        Maybe put wrapper in MLSFiles?
!          noSwaths = SWInqSwath ( fileNameString, allSwathNames, listSize )
          allSwathNames = ''
          noSwaths = mls_InqSwath ( fileNameString, allSwathNames, listSize )
          if ( listSize < len(allSwathNames) ) then
            commaPos = index ( allSwathNames, ',' )
            if ( commaPos == 0 ) commaPos = len_trim(allSwathNames)
            swathNameString = allSwathNames ( 1:commaPos )
          else
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Failed to determine swath names, string too long.' )
          end if
        endif
        
        ! Open the l2gp file
!        fileHandle = swopen(FileNameString, DFACC_READ)
        fileHandle = mls_io_gen_openF('swopen', .TRUE., returnStatus, &
             & record_length, DFACC_READ, FileName=FileNameString, &
             & debugOption=.true. )
        if ( fileHandle == -1 ) then
          call announce_error ( son, &
            & 'Failed to open swath file ' // trim(FileNameString) )
        end if

        ! Read the swath
        call ReadL2GPData ( fileHandle, swathNameString, l2gp )

        if( index(switches, 'apr') /= 0 ) &
        & call dump( l2gp, details=details )

        ! Close the file
!        fileHandle = swclose(fileHandle)
        fileHandle = mls_io_gen_closeF('swclose', fileHandle)
        if ( fileHandle == -1 ) then
          call announce_error ( son, &
            & 'Failed to close swath file ' // trim(FileNameString) )
        end if

        ! Add this l2gp to the database, decorate this key with index
        call decorate ( key, AddL2GPToDatabase( L2GPDatabase, l2gp ) )
        ! Don't call destroy contents as the AddL2GPToDatabase has done a shallow
        ! copy.
      case ( s_l2aux )

        if ( .not. all(got((/f_sdName, f_file/)))) &
          & call announce_error ( son, &
            & 'file/sd name must both be specified in read a priori' )
        
        call get_string ( sdName, sdNameString )
        sdNameString = sdNameString(2:LEN_TRIM(sdNameString)-1)
        ! create SD interface identifier for l2aux
        sd_id = sfstart(FilenameString, DFACC_READ)
        if (sd_id == -1 ) then
          call announce_error ( son, 'Failed to open l2aux ' // &
          &  trim(FilenameString) )
        end if

        ! Create the l2aux, and add it to the database.
        ! This doesn't match the interface in module L2AUXData
        !       CALL SetupNewL2AUXRecord ( l2aux )
        ! It has been relocated to READL2AUXData

        l2aux%name = l2Name

        l2Index = AddL2AUXToDatabase( L2AUXDatabase, l2aux )
        call decorate ( key, l2Index )
        call ReadL2AUXData ( sd_id, sdNameString, L2AUXDatabase(l2Index) )

        if( index(switches, 'apr') /= 0 ) &
        & call dump( L2AUXDatabase(l2Index), details )

      case ( s_gridded )

        if ( .not. all(got((/f_origin, f_field/))) ) &
          & call announce_error ( son, 'Incomplete gridded data information' )

        call get_string ( fieldName, fieldNameString, strip=.true. )
        
        select case ( griddedOrigin )
        case ( l_ncep )
          
          gridIndex = AddGriddedDataToDatabase( GriddedDatabase, GriddedData )
          call decorate ( key, gridIndex )
          call readGriddedData ( FileNameString, son, 'ncep', v_is_pressure, &
            & GriddedDatabase(gridIndex), &
            & 'XDim,YDim,Height,TIME', TRIM(fieldNameString) )
          
        case ( l_dao )
          
          gridIndex = AddGriddedDataToDatabase( GriddedDatabase, GriddedData )
          call decorate ( key, gridIndex )
          call ReadGriddedData ( FileNameString, son, 'dao', v_is_pressure, &
            & GriddedDatabase(gridIndex), fieldName = TRIM(fieldNameString) )
          
        case ( l_climatology )
          
          ! Identify file (maybe from PCF if no name given)

          if ( .NOT. got(f_file) .and. PCF) then
            
            do pcf_indx = lastClimPCF+1, mlspcf_l2clim_end
              returnStatus = Pgs_pc_getReference(i, version, fileNameString)
              if ( returnStatus == PGS_S_SUCCESS) exit
            end do
            
            if ( returnStatus /= PGS_S_SUCCESS ) then
              call announce_error ( son, &
                & 'PCF number not found to supply' // &
                & ' missing Climatology file name' )
            end if
          end if
          
          ! Have we read this already?
          gotAlready = associated(GriddedDatabase)
          if ( gotAlready ) then
            gotAlready = any(GriddedDatabase%sourceFilename==filenameString)
          end if
          if ( .not. gotAlready ) then
            ! No, well read it then, add its entire contents to the database
            call read_climatology ( FileNameString, son, &
              & GriddedDatabase, mlspcf_l2clim_start, mlspcf_l2clim_end )
          end if
       
          ! Locate requested grid by name, store index in gridIndex
          
          ! Check that field name is among those added by the source field
          do gridIndex = 1, size(griddedDatabase)
            if ( trim(fieldNameString) == &
              & trim(GriddedDatabase(gridIndex)%quantityName) ) exit
          end do

          if ( gridIndex <= size(griddedDatabase) ) then
            call decorate ( key, gridIndex )
          else
            call announce_error ( son, 'Field ' // trim(fieldNameString) // &
              & ' not found in clim. file ' // trim(fileNameString) )
          end if
        case default ! Can't get here if tree_checker worked correctly
        end select   ! origins of gridded data

        if( index(switches, 'apr') /= 0 ) &
        & call dump( GriddedDatabase(gridIndex), details )

      case default
      end select     ! types of apriori data
      
    end do                              ! Lines in l2cf loop
    
    if ( ERROR/=0 ) then
      call MLSMessage(MLSMSG_Error,ModuleName, &
        & 'Problem with read_apriori section')
    end if

    if ( toggle(gen) ) call trace_end("read_apriori")
  
    if ( timing ) call sayTime
    return

  contains
    subroutine SayTime
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      endif
      call output ( "Timing for read_apriori = " )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime
  end subroutine read_apriori

  ! ------------------------------------------------  announce_error  -----
  subroutine Announce_error ( lcf_where, full_message, use_toolkit, &
    & error_number )
  
   ! Arguments
  
    integer, intent(in)    :: Lcf_where
    character(LEN=*), intent(in)    :: Full_message
    logical, intent(in), optional :: Use_toolkit
    integer, intent(in), optional    :: Error_number
    ! Local
    logical :: Just_print_it
    logical, parameter :: Default_output_by_toolkit = .true.
 
    if ( present(use_toolkit) ) then
      just_print_it = .not. use_toolkit
    else if ( default_output_by_toolkit ) then
      just_print_it = .false.
    else
      just_print_it = .true.
    end if
 
    if ( .not. just_print_it ) then
      error = max(error,1)
      call output ( '***** At ' )

      if ( lcf_where > 0 ) then
          call print_source ( source_ref(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      end if

      call output ( ': ' )
      call output ( "The " );
      if ( lcf_where > 0 ) then
        call dump_tree_node ( lcf_where, 0 )
      else
        call output ( '(no lcf tree available)' )
      end if

      call output ( " Caused the following error: ", advance='yes', &
        & from_where=ModuleName )
      call output ( trim(full_message), advance='yes', &
        & from_where=ModuleName )
      if ( present(error_number) ) then
        call output ( 'error number ', advance='no' )
        call output ( error_number, places=9, advance='yes' )
      end if
    else
      call output ( '***Error in module ' )
      call output ( ModuleName, advance='yes' )
      call output ( trim(full_message), advance='yes' )
      if ( present(error_number) ) then
        call output ( 'Error number ' )
        call output ( error_number, advance='yes' )
      end if
    end if

!===========================
  end subroutine Announce_error
!===========================

end module ReadAPriori

!=============================================================================

!
! $Log$
! Revision 2.28  2002/01/18 00:55:30  pwagner
! Uses MLSFiles for swapi wrappers
!
! Revision 2.27  2002/01/09 00:00:04  pwagner
! Replaced write or print statements with calls to output
!
! Revision 2.26  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.25  2001/10/30 00:35:27  pwagner
! Tidied up small things
!
! Revision 2.24  2001/10/26 23:20:10  pwagner
! Optionally dumps apriori quantities as it reads them
!
! Revision 2.23  2001/10/08 21:35:43  pwagner
! Initialize allswathnames
!
! Revision 2.22  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.21  2001/09/10 23:37:44  livesey
! New GriddedData stuff
!
! Revision 2.20  2001/05/12 00:20:00  livesey
! Allowed user to not supply swath name when reading l2gp
!
! Revision 2.19  2001/05/07 18:03:56  pwagner
! Checks for PCF before looping over pcf_indx
!
! Revision 2.18  2001/05/03 20:34:08  vsnyder
! Cosmetic changes
!
! Revision 2.17  2001/04/16 23:50:01  pwagner
! Tiny change to announce_error
!
! Revision 2.16  2001/04/12 22:19:33  vsnyder
! Improved an error message
!
! Revision 2.15  2001/04/10 20:04:26  livesey
! Bug fixes etc.
!
! Revision 2.14  2001/03/30 00:27:38  pwagner
! Fleshed out njl outline for reading clim. files
!
! Revision 2.13  2001/03/29 19:13:41  livesey
! Added some comment guidelines for read climatology
!
! Revision 2.12  2001/03/21 00:46:08  pwagner
! Passes son to READ_CLIMATOLOGY
!
! Revision 2.11  2001/03/15 21:38:00  pwagner
! Gets v_is_pressure from GriddedData
!
! Revision 2.10  2001/03/15 21:25:16  pwagner
! Split between GriddedData and ncep_dao modules
!
! Revision 2.9  2001/03/15 21:18:57  vsnyder
! Use Get_Spec_ID instead of decoration(subtree...
!
! Revision 2.8  2001/03/15 00:34:34  pwagner
! Gives more info to ReadGriddedData
!
! Revision 2.7  2001/03/14 19:06:08  livesey
! Some changes
!
! Revision 2.6  2001/03/14 18:54:38  pwagner
! Uses FieldNameString and son in call to ReadGriddedData
!
! Revision 2.5  2001/03/08 01:08:08  pwagner
! Interfaces with ReadGriddedData
!
! Revision 2.4  2001/03/07 22:41:44  livesey
! Reworked the l2cf aspects
!
! Revision 2.3  2001/03/07 01:04:33  pwagner
! No longer uses obtainclim, obtaindao, obtainncep
!
! Revision 2.2  2001/03/06 00:23:58  pwagner
! A little bit more
!
! Revision 2.1  2001/03/03 00:14:40  pwagner
! First commit
!
!
