! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ReadAPriori

  use GriddedData, only: GriddedData_T, AddGridTemplateToDatabase, &
  & READ_CLIMATOLOGY, ReadGriddedData
  use Hdf, only: DFACC_READ, SFSTART
  use Hdfeos, only: swopen, swclose
  use INIT_TABLES_MODULE, only: F_FIELD, F_FILE, F_ORIGIN, F_SDNAME, F_SWATH, &
    & L_CLIMATOLOGY, L_DAO, L_NCEP, S_GRIDDED, S_L2AUX, S_L2GP, LIT_INDICES
  use L2AUXData, only: L2AUXData_T, AddL2AUXToDatabase, &
    &                  ReadL2AUXData!, SetupNewL2AUXRecord
  use L2GPData, only: L2GPData_T, AddL2GPToDatabase, ReadL2GPData
  use MLSCommon, only: FileNameLen, L1BInfo_T, TAI93_Range_T
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    &                         MLSMSG_Error, MLSMSG_FileOpen!, MLSMSG_Info
!  use OBTAINCLIMATOLOGY, only: READ_CLIMATOLOGY
!  use OBTAINDAO, only: READ_DAO
!  use OBTAINNCEP, only: READ_NCEP
  use MLSPCF2, only: mlspcf_l2clim_start, mlspcf_l2clim_end
  use String_Table, only: GET_STRING, DISPLAY_STRING
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, &
    &             SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_NAMED!, N_DOT
  use OUTPUT_M, only: OUTPUT

  implicit none
  private
  public ::  read_apriori

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

  subroutine read_apriori ( root, L2GPDatabase, l2auxDatabase, aprioriData)

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the Read a priori section in the AST
    type (l2gpdata_t), dimension(:), pointer :: L2GPDatabase
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase
    type (GriddedData_T), dimension(:), pointer :: aprioriData 

    !Local Variables
    integer :: FIELD               ! Son of KEY, must be n_assign
    integer :: fileHandle          ! fileHandle of a priori data file
    integer :: fileName            ! Sub-rosa index of name in file='name'
    character(len=FileNameLen) :: FileNameString   ! actual literal file name
    integer :: FileType            ! either s_l2gp or s_l2aux
    integer :: I, J                ! Loop indices for section, spec
    integer :: KEY                 ! Index of n_spec_args in the AST
    type (L2AUXData_T) :: L2AUX
    type (L2GPData_T) :: L2GP
    TYPE (GriddedData_T) :: GriddedData
    integer :: l2Index             ! In the l2gp or l2aux database
    integer :: L2Name              ! Sub-rosa index of L2[aux/gp] label
    character (LEN=480) :: msr     ! Error message if can't find file

    integer :: sd_id
    integer :: SON              ! Of root, an n_spec_args or a n_named
    integer :: swathName        ! sub-rosa index of name in swath='name'
    integer :: sdName        ! sub-rosa index of name in sdName='name'
    integer :: fieldName        ! sub-rosa index of name in field='name'
    integer :: griddedOrigin            ! From tree
    character(len=FileNameLen) :: SWATHNAMESTRING ! actual literal swath name
    character(len=FileNameLen) :: SDNAMESTRING ! actual literal sdName
    character(len=FileNameLen) :: FIELDNAMESTRING ! actual literal fieldName

    if ( toggle (gen) ) call trace_begin( "read_apriori", root )

    do i = 2, nsons(root)-1 ! Skip the section name at begin and end
      son = subtree(i,root)
      if ( node_id(son) == n_named ) then ! Is spec labeled?
        key = subtree(2,son)
        l2Name = sub_rosa(subtree(1,son))
      else
        key = son
        l2Name = 0
      end if

      ! Node_id(key) is now n_spec_args.

      FileType = decoration(subtree(1,decoration(subtree(1,key))))

      ! Now parse file and field names
      fileName = 0
      swathName = 0
      do j = 2, nsons(key)
        field = subtree(j,key)
        select case ( decoration(subtree(1,field)) )
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

      if ( fileName == 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'File name not specified in read a priori')
      
      call get_string ( FileName, fileNameString )
      fileNameString=fileNameString(2:LEN_TRIM(fileNameString)-1)
      
      select case( FileType )
      case ( s_l2gp )

        if ( swathName == 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
          & 'Swath name not specified in read a priori')
        
        call get_string ( swathName, swathNameString )
        swathNameString=swathNameString(2:LEN_TRIM(swathNameString)-1)
        ! Open the l2gp file
        fileHandle = swopen(FileNameString, DFACC_READ)
        if (fileHandle == -1) then
          msr = MLSMSG_Fileopen // FileNameString
          call MLSMessage ( MLSMSG_Error, ModuleName, trim(msr) )
        end if

        ! Read the swath
        call ReadL2GPData ( fileHandle, swathNameString, l2gp )

        ! Close the file
        fileHandle = swclose(fileHandle)
        if (fileHandle == -1) THEN
          msr = 'Failed to close file ' // FileNameString
          call MLSMessage(MLSMSG_Error, ModuleName, trim(msr))
        end if

        ! Add this l2gp to the database, decorate this key with index
        call decorate ( key, AddL2GPToDatabase( L2GPDatabase, l2gp ) )
        ! Don't call destroy contents as the AddL2GPToDatabase has done a shallow
        ! copy.
      case ( s_l2aux )

        if ( sdName == 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
          & 'sd name not specified in read a priori')
        
        call get_string ( sdName, sdNameString )
        sdNameString=sdNameString(2:LEN_TRIM(sdNameString)-1)
        ! create SD interface identifier for l2aux
        sd_id = sfstart(FilenameString, DFACC_READ)
        IF (sd_id == -1) THEN
          msr = MLSMSG_Fileopen // FileNameString
          CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
        ENDIF
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
        CALL ReadL2AUXData(sd_id, sdNameString, L2AUXDatabase(l2Index))

      case ( s_gridded )

        call output('Hello paul, in read gridded data section.',advance='yes')
        call display_string(lit_indices(griddedOrigin),advance='yes')
        if ( fieldName == 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
          & 'Field name not specified in read a priori')

        call get_string ( fieldName, fieldNameString )
        fieldNameString=fieldNameString(2:LEN_TRIM(fieldNameString)-1)
        
        select case ( griddedOrigin )
        case ( l_ncep )
          
          l2Index = AddGridTemplateToDatabase( aprioriData, GriddedData )
          call decorate ( key, l2Index )
          CALL ReadGriddedData(FileNameString, &
            & aprioriData(l2Index), 'Some_field_name')
          
        case ( l_dao )
          
          l2Index = AddGridTemplateToDatabase( aprioriData, GriddedData )
          call decorate ( key, l2Index )
          CALL ReadGriddedData(FileNameString, &
            & aprioriData(l2Index), 'Some_field_name')
          
        case ( l_climatology )
          
          CALL READ_CLIMATOLOGY(FileNameString, &
            & aprioriData, mlspcf_l2clim_start, mlspcf_l2clim_end)
          
        case default ! Can't get here if tree_checker worked correctly
        end select
      case default
      end select
      
    end do                              ! Lines in l2cf loop
    
    if (toggle(gen) ) call trace_end("read_apriori")
  
  end subroutine read_apriori

end module ReadAPriori

!=============================================================================

!
! $Log$
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
