! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ReadAPriori

  use GriddedData, only: GriddedData_T, v_is_pressure
  use Hdf, only: DFACC_READ, SFSTART
  use Hdfeos, only: swopen, swclose
  use INIT_TABLES_MODULE, only: F_FIELD, F_FILE, F_ORIGIN, F_SDNAME, F_SWATH, &
    & FIELD_FIRST, FIELD_LAST, L_CLIMATOLOGY, L_DAO, L_NCEP, S_GRIDDED, S_L2AUX, &
    & S_L2GP, LIT_INDICES
  use L2AUXData, only: L2AUXData_T, AddL2AUXToDatabase, &
    &                  ReadL2AUXData!, SetupNewL2AUXRecord
  use L2GPData, only: L2GPData_T, AddL2GPToDatabase, ReadL2GPData
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: FileNameLen, L1BInfo_T, TAI93_Range_T
  use MLSFiles, only: SPLIT_PATH_NAME
!  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
!    &                         MLSMSG_Error, MLSMSG_FileOpen!, MLSMSG_Info
!  use OBTAINCLIMATOLOGY, only: READ_CLIMATOLOGY
!  use OBTAINDAO, only: READ_DAO
!  use OBTAINNCEP, only: READ_NCEP
  use MLSPCF2, only: mlspcf_l2clim_start, mlspcf_l2clim_end
  use MoreTree, only: Get_Spec_ID
  use ncep_dao, only: AddGridTemplateToDatabase, &
  & READ_CLIMATOLOGY, ReadGriddedData, source_file_already_read
  use OUTPUT_M, only: OUTPUT
  use SDPToolkit, only: Pgs_pc_getReference, PGS_S_SUCCESS
  use String_Table, only: GET_STRING, DISPLAY_STRING
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, &
    &             SUB_ROSA, SUBTREE, DUMP_TREE_NODE, SOURCE_REF
  use TREE_TYPES, only: N_NAMED!, N_DOT

  implicit none
  private
  public ::  read_apriori
  private ::  announce_error
  integer, private :: ERROR

  ! -----     Private declarations     ---------------------------------

  private :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  character(len=130) :: id = &
    "$id: open_init.f90,v 1.11 2000/06/19 22:40:51 lungu Exp $"
  character(len=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! --------------------------------------------------  read_apriori  -----
  ! Read a priori data from data files, be they l2gp, l2aux, climatology,
  ! NCEP, DAO etc.

  subroutine read_apriori ( root, L2GPDatabase, l2auxDatabase, GriddedDatabase)

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the Read a priori section in the AST
    type (l2gpdata_t), dimension(:), pointer :: L2GPDatabase
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase
    type (GriddedData_T), dimension(:), pointer :: GriddedDatabase 

    !Local Variables
    integer :: FIELD               ! Son of KEY, must be n_assign
    integer :: FIELDINDEX          ! Literal
    integer :: fileHandle          ! fileHandle of a priori data file
    integer :: fileName            ! Sub-rosa index of name in file='name'
    character(len=FileNameLen) :: FileNameString   ! actual literal file name
    integer :: FileType            ! either s_l2gp or s_l2aux
    integer :: gridIndex           ! In the griddeddata database
    logical :: gotAlready               ! Do we need to reread this file?
    integer :: I, J                ! Loop indices for section, spec
    integer :: KEY                 ! Index of n_spec_args in the AST
    integer :: lastClimPCF
    type (L2AUXData_T) :: L2AUX
    type (L2GPData_T) :: L2GP
    TYPE (GriddedData_T) :: GriddedData
    integer :: l2Index             ! In the l2gp or l2aux database
    integer :: L2Name              ! Sub-rosa index of L2[aux/gp] label
    character (LEN=480) :: msr     ! Error message if can't find file
    integer :: pcf                 ! loop index of climatology pcf numbers

    integer :: sd_id
    integer :: SON              ! Of root, an n_spec_args or a n_named
    integer :: swathName        ! sub-rosa index of name in swath='name'
    integer :: sdName        ! sub-rosa index of name in sdName='name'
    integer :: fieldName        ! sub-rosa index of name in field='name'
    integer :: griddedOrigin            ! From tree
    character(len=FileNameLen) :: FIELDNAMESTRING ! actual literal clim. field
    character(len=FileNameLen) :: SWATHNAMESTRING ! actual literal swath name
    character(len=FileNameLen) :: SDNAMESTRING ! actual literal sdName
    character(len=FileNameLen) :: bareFilename
    character(len=FileNameLen) :: path
    integer :: version
	integer :: returnStatus

    logical, dimension(field_first:field_last) :: got

    if ( toggle (gen) ) call trace_begin( "read_apriori", root )

    error = 0
	 version = 1
	lastClimPCF = mlspcf_l2clim_start - 1

    do i = 2, nsons(root)-1 ! Skip the section name at begin and end
      got=.false.
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
        fieldIndex=decoration(subtree(1,field))
        got(fieldIndex)=.true.
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

!      if ( fileName == 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
!        & 'File name not specified in read a priori')
      
      if ( got(f_file) ) then
        call get_string ( FileName, fileNameString, strip=.true. )
      else
        fileNameString=''
      end if
        
      select case( FileType )
      case ( s_l2gp )
!        if ( .not. all(got((/f_swath, f_file/)))) &
!          & call MLSMessage(MLSMSG_Error, ModuleName, &
!          & 'Swath/filename name not specified in read a priori')
        if ( .not. all(got((/f_swath, f_file/)))) &
          & call announce_error ( son, &
            & 'Swath/filename name not specified in read a priori' )
        
        call get_string ( swathName, swathNameString )
        swathNameString = swathNameString(2:LEN_TRIM(swathNameString)-1)
        ! Open the l2gp file
        fileHandle = swopen(FileNameString, DFACC_READ)
        if ( fileHandle == -1 ) then
!          msr = MLSMSG_Fileopen // FileNameString
!          call MLSMessage ( MLSMSG_Error, ModuleName, trim(msr) )
          call announce_error ( son, &
            & 'Failed to open swath file ' // FileNameString )
        end if

        ! Read the swath
        call ReadL2GPData ( fileHandle, swathNameString, l2gp )

        ! Close the file
        fileHandle = swclose(fileHandle)
        if ( fileHandle == -1 ) then
!          msr = 'Failed to close file ' // FileNameString
!          call MLSMessage(MLSMSG_Error, ModuleName, trim(msr))
          call announce_error ( son, &
            & 'Failed to close swath file ' // FileNameString )
        end if

        ! Add this l2gp to the database, decorate this key with index
        call decorate ( key, AddL2GPToDatabase( L2GPDatabase, l2gp ) )
        ! Don't call destroy contents as the AddL2GPToDatabase has done a shallow
        ! copy.
      case ( s_l2aux )

!        if ( .not. all(got((/f_sdName, f_file/)))) &
!          & call MLSMessage(MLSMSG_Error, ModuleName, &
!          & 'file/sd name not specified in read a priori')
        if ( .not. all(got((/f_sdName, f_file/)))) &
          & call announce_error ( son, &
            & 'file/sd name not specified in read a priori' )
        
        call get_string ( sdName, sdNameString )
        sdNameString = sdNameString(2:LEN_TRIM(sdNameString)-1)
        ! create SD interface identifier for l2aux
        sd_id = sfstart(FilenameString, DFACC_READ)
        if (sd_id == -1 ) then
!          msr = MLSMSG_Fileopen // FileNameString
!          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
          call announce_error ( son, 'Failed to open l2aux ' // sdNameString )
        end if
        ! ??? subtree(1,key) is l2aux or l2gp.  It doesn't have a subtree ???
        !       vectorIndex = decoration(decoration(subtree(2,subtree(1,key))))

        ! Create the l2aux, and add it to the database.
        ! This doesn't match the interface in module L2AUXData
        !       CALL SetupNewL2AUXRecord ( l2aux )
        ! It has been relocated to READL2AUXData

        l2aux%name = l2Name

        !        call decorate ( key, AddL2AUXToDatabase( L2AUXDatabase, l2aux ) )

        ! That's the end of the create operation

        ! ??? Should "vectorIndex" be "decoration(key)" ???
        ! ??? If so, do something like
        l2Index = AddL2AUXToDatabase( L2AUXDatabase, l2aux )
        call decorate ( key, l2Index )
        !   call ReadL2AUXData ( ... L2AUXDataBase(l2Index) ... )
        ! Need to add this routine to L2AUXData.f90 before uncommenting this line
        call ReadL2AUXData(sd_id, sdNameString, L2AUXDatabase(l2Index))

      case ( s_gridded )

!        if ( .not. all(got((/f_origin, f_field, f_file/)))) &
!          & call MLSMessage(MLSMSG_Error, ModuleName, &
!          & 'Incomplete gridded data information')
        if ( .not. all(got((/f_origin, f_field/)))) &
          & call announce_error(son, &
          & 'Incomplete gridded data information')

        call get_string ( fieldName, fieldNameString, strip=.true. )
        
        select case ( griddedOrigin )
        case ( l_ncep )
          
          gridIndex = AddGridTemplateToDatabase( GriddedDatabase, GriddedData )
          call decorate ( key, gridIndex )
          CALL ReadGriddedData ( FileNameString, son, 'ncep', v_is_pressure, &
            & GriddedDatabase(gridIndex), &
            & 'XDim,YDim,Height,TIME', TRIM(fieldNameString) )
          
        case ( l_dao )
          
          gridIndex = AddGridTemplateToDatabase( GriddedDatabase, GriddedData )
          call decorate ( key, gridIndex )
          call ReadGriddedData ( FileNameString, son, 'dao', v_is_pressure, &
            & GriddedDatabase(gridIndex), fieldName = TRIM(fieldNameString) )
          
        case ( l_climatology )
          
          ! Identify file (maybe from PCF if no name given)

          if ( .NOT. got(f_file) ) then
            
            do pcf=lastClimPCF+1, mlspcf_l2clim_end
              returnStatus = Pgs_pc_getReference(i, version, &
                & fileNameString)
              if ( returnStatus == PGS_S_SUCCESS) exit
            end do
            
            if ( returnStatus /= PGS_S_SUCCESS ) then
              call announce_error ( son, &
                & 'PCF number not found to supply missing Climatology file name' )
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
          gridIndex=1
          do gridIndex=1,size(griddedDatabase)
            if ( trim(fieldNameString) == &
              & trim(GriddedDatabase(gridIndex)%quantityName)) exit
          end do

          if ( gridIndex <= size(griddedDatabase) ) then
            call decorate ( key, gridIndex )
          else
            call announce_error ( son, 'Field ' // trim(fieldNameString) // &
              & ' not found in clim. file ' // trim(fileNameString) )
          end if
        case default ! Can't get here if tree_checker worked correctly
        end select
      case default
      end select
      
    end do                              ! Lines in l2cf loop
    
    if ( toggle(gen) ) call trace_end("read_apriori")
  
  end subroutine read_apriori

  ! ------------------------------------------------  announce_error  -----
  subroutine announce_error ( lcf_where, full_message, use_toolkit, &
  & error_number )
  
   ! Arguments
  
    integer, intent(in)    :: lcf_where
    character(LEN=*), intent(in)    :: full_message
    logical, intent(in), optional :: use_toolkit
    integer, intent(in), optional    :: error_number
    ! Local
!    character (len=80) :: msg, mnemonic
!    integer :: status
    logical :: just_print_it
    logical, parameter :: default_output_by_toolkit = .true.
 
    if ( present(use_toolkit) ) then
      just_print_it = use_toolkit
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
      print*, '***Error in module ', ModuleName
      print*, trim(full_message)
      if ( present(error_number) ) then
        print*, 'error number ', error_number
      end if
    end if

!===========================
  end subroutine announce_error
!===========================

end module ReadAPriori

!=============================================================================

!
! $Log$
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
