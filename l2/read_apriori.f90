! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ReadAPriori

  use GriddedData, only: GriddedData_T, AddGridTemplateToDatabase
  use Hdf, only: DFACC_READ, SFSTART
  use Hdfeos, only: swopen, swclose
  use INIT_TABLES_MODULE, only: F_FILE, F_SWATH, S_L2AUX, S_L2GP !, s_ncep, &
	!  & s_dao, s_clim 
  use L2AUXData, only: L2AUXData_T, AddL2AUXToDatabase, &
    &                  ReadL2AUXData!, SetupNewL2AUXRecord
  use L2GPData, only: L2GPData_T, AddL2GPToDatabase, ReadL2GPData
  use MLSCommon, only: FileNameLen, L1BInfo_T, TAI93_Range_T
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    &                         MLSMSG_Error, MLSMSG_FileOpen!, MLSMSG_Info
  use OBTAINCLIMATOLOGY, only: READ_CLIMATOLOGY
  use OBTAINDAO, only: READ_DAO
  use OBTAINNCEP, only: READ_NCEP
  use String_Table, only: GET_STRING
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, &
    &             SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_NAMED!, N_DOT

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
    integer :: l2Index             ! In the l2gp or l2aux database
    integer :: L2Name              ! Sub-rosa index of L2[aux/gp] label
    character (LEN=480) :: msr     ! Error message if can't find file
!?  type (Vector_T) :: newVector

    integer :: sd_id
    integer :: SON              ! Of root, an n_spec_args or a n_named
    integer :: swathName        ! sub-rosa index of name in swath='name'
    character(len=FileNameLen) :: SwathNameString ! actual literal swath name
!?  integer :: vectorIndex         ! In the vector database


! These should be replaced by appropriate entries in init_tables_module
	INTEGER, PARAMETER :: s_ncep=s_l2aux+1
	INTEGER, PARAMETER :: s_dao=s_ncep+1
	INTEGER, PARAMETER :: s_climatology=s_dao+1

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
        end select
      end do
      if ( fileName == 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'File name not specified in read a priori')
      if ( swathName == 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Swath name not specified in read a priori')

      call get_string ( FileName, fileNameString )
      call get_string ( swathName, swathNameString )
      fileNameString=fileNameString(2:LEN_TRIM(fileNameString)-1)
      swathNameString=swathNameString(2:LEN_TRIM(swathNameString)-1)

      select case( FileType )
      case ( s_l2gp )

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
        CALL ReadL2AUXData(sd_id, swathNameString, L2AUXDatabase(l2Index))

	! The remaining cases are gridded data types
      case ( s_ncep )

			CALL READ_NCEP(FileNameString, &
			& aprioriData(l2Index)%field(1, 1, 1, :, :, :))

      case ( s_dao )

			CALL READ_DAO(FileNameString, 'Some_vector_name', &
			& aprioriData(l2Index)%field(1, 1, 1, :, :, :))

      case ( s_climatology )

			CALL READ_CLIMATOLOGY(FileNameString, &
			& aprioriData(l2Index)%field(1, 1, 1, :, :, :))

      case default ! Can't get here if tree_checker worked correctly
      end select


    end do

    if (toggle(gen) ) call trace_end("read_apriori")

  end subroutine read_apriori

end module ReadAPriori

!=============================================================================

!
! $Log$
! Revision 2.1  2001/03/03 00:14:40  pwagner
! First commit
!
!
