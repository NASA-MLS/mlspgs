! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Open_Init

  ! Opens and closes several files
  ! Creates and destroys the L1BInfo database

  use Hdf, only: DFACC_READ, SFSTART
  use Hdfeos, only: swopen, swclose
  use INIT_TABLES_MODULE, only: F_FILE, F_SWATH, S_L2AUX, S_L2GP
  use L2AUXData, only: AddL2AUXToDatabase, L2AUXData_T, ReadL2AUXData
  use L2GPData, only: AddL2GPToDatabase, L2GPData_T, ReadL2GPData
  use MLSCommon, only: FileNameLen, L1BInfo_T, TAI93_Range_T
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    &                         MLSMSG_Error, MLSMSG_FileOpen!, MLSMSG_Info
  use MLSPCF2, only: MLSPCF_L1B_OA_START, MLSPCF_L1B_RAD_END, &
    &                MLSPCF_L1B_RAD_START, MLSPCF_NOMEN_START!, &
!   &                MLSPCF_L2CF_START
  use MLSSignalNomenclature, only: ReadSignalsDatabase
  use SDPToolkit, only: PGS_IO_Gen_closeF, PGS_IO_Gen_openF, &
    &                   Pgs_pc_getReference, PGS_S_SUCCESS, &
    &                   PGSd_IO_Gen_RSeqFrm, PGSTD_E_NO_LEAP_SECS
  use String_Table, only: Get_String !, L2CFUnit => INUNIT
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
! use TREE, only: DECORATE, DECORATION, DUMP_TREE_NODE, NODE_ID, NSONS, &
!   &             SOURCE_REF, SUB_ROSA, SUBTREE
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, &
    &             SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_NAMED!, N_DOT
! use VectorsModule, only: AddVectorToDatabase, CreateVector, Dump, Vector_T, &
!   &                      VectorTemplate_T
! use VectorsModule, only: Vector_T

  implicit none
  private
  public :: Close_MLSCF, DestroyL1BInfo, OpenAndInitialize, Open_MLSCF
  public :: Read_apriori

  ! -----     Private declarations     ---------------------------------

  private :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  character(len=130) :: id = &
    "$id: open_init.f90,v 1.11 2000/06/19 22:40:51 lungu Exp $"
  character(len=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! ------------------------------------------------  Close_MLSCF  -----
  subroutine Close_MLSCF ( CF_Unit )

    integer, intent(in) :: CF_Unit

    character (len=32) :: Mnemonic
    character (len=256) :: Msg
    integer :: ReturnStatus

    returnStatus = Pgs_io_gen_closeF ( CF_Unit )

    if ( returnStatus /= PGS_S_SUCCESS ) then
      call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Error closing L2CF:  '//mnemonic//' '//msg)
    end if

  end subroutine Close_MLSCF

  ! ---------------------------------------------  DestroyL1BInfo  -----
  subroutine DestroyL1BInfo ( L1BInfo )
    type (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset
    integer :: STATUS ! from deallocate
    if (associated(l1bInfo%L1BRADIDs)) then
       deallocate( l1bInfo%L1BRADIDs, stat=status )
       if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & MLSMSG_DeAllocate // "l1bInfo" )
    endif
  end subroutine DestroyL1BInfo


  ! ------------------------------------------  OpenAndInitialize  -----
  subroutine OpenAndInitialize ( processingRange, l1bInfo )

    ! Opens L1 RAD files
    ! Opens L1OA file
    ! Opens and reads the signal nomenclature file
    ! Gets the start and end times from the PCF

    ! Arguments

    type (TAI93_Range_T) :: processingRange ! Data processing range
    type (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset

    !Local Variables
    integer, parameter :: CCSDSEndId = 10412
    integer, parameter :: CCSDSLen=27
    integer, parameter :: CCSDSStartId = 10411

    character(len=CCSDSlen) CCSDSEndTime
    character(len=CCSDSlen) CCSDSStartTime
    integer :: ifl1
    integer :: L1FileHandle, L1_Version
    character (LEN=132) :: L1physicalFilename
    integer :: nomenUnit
    integer :: Nomen_Version
    integer :: returnStatus
    integer :: sd_id
    integer :: STATUS ! From allocate

    integer :: pgs_td_utctotai, pgs_pc_getconfigdata

    ifl1 = 0

    ! Open L1 RAD files
    ! Get the l1 file name from the PCF
    L1_Version = 1

    do L1FileHandle = mlspcf_l1b_rad_start, mlspcf_l1b_rad_end

      returnStatus = Pgs_pc_getReference(L1FileHandle, L1_Version, &
        & L1physicalFilename)

      if (returnStatus == PGS_S_SUCCESS) then

        ! Open the HDF file and initialize the SD interface

        ! Allocate L1BRADIDs

        allocate( l1bInfo%L1BRADIDs(10), stat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Allocate // "l1bInfo" )

        sd_id = sfstart(L1physicalFilename, DFACC_READ)
        if (sd_id == -1) then
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Error opening L1RAD file "//L1physicalFilename )

        else
          ifl1 = ifl1 + 1
          l1bInfo%L1BRADIDs(ifl1) = sd_id
        end if
      end if
    end do ! L1FileHandle = mlspcf_l1b_rad_start, mlspcf_l1b_rad_end

    !if (ifl1 == 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
    !  & "Could not find any L1BRAD files" )

    ! Open L1OA File

    L1_Version = 1
    returnStatus = Pgs_pc_getReference(mlspcf_l1b_oa_start, L1_Version, &
      & L1physicalFilename)

    if (returnStatus == PGS_S_SUCCESS) then

      ! Open the HDF file and initialize the SD interface

      sd_id = sfstart(L1physicalFilename, DFACC_READ)
      if (sd_id == -1) then

        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Error opening L1OA file "//L1physicalFilename )

      else
        l1bInfo%L1BOAID = sd_id
      end if

    else

      call MLSMessage ( MLSMSG_Error, ModuleName, "Could not find L1BOA file" )

    end if

    ! Open and Read nomenclature file

    Nomen_Version = 1
    returnStatus = PGS_IO_Gen_openF (mlspcf_nomen_start, PGSd_IO_Gen_RSeqFrm, &
      & 0, NomenUnit, Nomen_Version)
    if ( returnStatus /= PGS_S_SUCCESS ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "Could not open nomenclature file" )

    call ReadSignalsDatabase ( Nomenunit )

    returnStatus = PGS_IO_Gen_closeF (Nomenunit)  ! close unit
    if ( returnstatus /= PGS_S_SUCCESS )  call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "Could not close nomenclature file" )

    ! Get the Start and End Times from PCF

    returnStatus = pgs_pc_getconfigdata (CCSDSStartId, CCSDSStartTime)
    if ( returnstatus /= PGS_S_SUCCESS ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "Could not get CCSDS Start Time" )

    returnStatus = pgs_td_utctotai (CCSDSStartTime, processingrange%starttime)
    !   ??? Is PGSTD_E_NO_LEAP_SECS an OK status ???
    if ( returnstatus /= PGS_S_SUCCESS .and. &
      &  returnstatus /= PGSTD_E_NO_LEAP_SECS ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "Could not convert UTC Start time to TAI" )

    returnStatus = pgs_pc_getconfigdata (CCSDSEndId, CCSDSEndTime)
    if ( returnstatus /= PGS_S_SUCCESS ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "Could not get CCSDS End Time" )

    returnStatus = pgs_td_utctotai (CCSDSEndTime, processingrange%endtime)
    !   ??? Is PGSTD_E_NO_LEAP_SECS an OK status ???
    if ( returnstatus /= PGS_S_SUCCESS .and. &
      & returnstatus /= PGSTD_E_NO_LEAP_SECS) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "Could not convert UTC Start time to TAI" )

    return

  end subroutine OpenAndInitialize

  ! -------------------------------------------------  Open_MLSCF  -----
  subroutine Open_MLSCF ( MLSPCF_Start, CF_Unit )

    integer, intent(in) :: MLSPCF_Start
    integer, intent(out) :: CF_Unit

    integer :: L2CF_Version
    character (len=32) :: Mnemonic
    character (len=256) :: Msg
    integer :: ReturnStatus

    !   Open the MLSCF as a generic file for reading
    L2CF_Version = 1
    returnStatus = Pgs_io_gen_openF ( mlspcf_start, PGSd_IO_Gen_RSeqFrm, &
      &                               0, CF_Unit, L2CF_Version)

    if ( returnStatus /= PGS_S_SUCCESS ) then
      call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "Error opening MLSCF:  "//mnemonic//"  "//msg)
    end if

  end subroutine Open_MLSCF

  ! --------------------------------------------------  read_apriori  -----
  ! Read a priori data from data files, be they l2gp, l2aux, climatology,
  ! NCEP, DAO etc.

  subroutine read_apriori ( root, L2GPDatabase, l2auxDatabase)


    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the Read a priori section in the AST
    type (l2gpdata_t), dimension(:), pointer :: L2GPDatabase
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase

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

      case default ! Can't get here if tree_checker worked correctly
      end select


    end do

    if (toggle(gen) ) call trace_end("read_apriori")

  end subroutine read_apriori

end module Open_Init
!=============================================================================

!
! $Log$
! Revision 2.21  2001/03/10 07:07:58  livesey
! Made it not mind if no L1B radiance files.
!
! Revision 2.20  2001/03/07 22:49:17  vsnyder
! Commented-out more USEd entities that NAG says actually aren't used.
!
! Revision 2.19  2001/03/03 00:08:58  pwagner
! Lost read_apriori and read_mlscf to new modules
!
! Revision 2.18  2001/02/27 01:30:46  vsnyder
! Commented-out several USEd entities that NAG says actually aren't used.
!
! Revision 2.17  2001/02/23 18:17:35  livesey
! Added trace calls
!
! Revision 2.16  2001/02/23 00:53:07  vsnyder
! Correct an error message
!
! Revision 2.15  2001/02/16 00:48:07  livesey
! Added stuff to read l2gp's in
!
! Revision 2.14  2001/02/13 22:59:36  pwagner
! l2 modules can only use MLSPCF2
!
! Revision 2.13  2001/02/09 00:38:22  livesey
! Various updates
!
! Revision 2.12  2001/02/08 00:58:14  vsnyder
! Correct calculation of "field"
!
! Revision 2.11  2001/01/03 00:47:21  pwagner
! Calls READL2AUXData from L2AUXData module
!
! Revision 2.10  2000/12/04 21:47:46  pwagner
! Uses parser better
!
! Revision 2.9  2000/12/02 01:11:59  pwagner
! Added ReadL2AUXData
!
! Revision 2.8  2000/12/02 00:00:40  vsnyder
! More misc. cleanup.
!
! Revision 2.7  2000/12/01 23:35:25  vsnyder
! Use abstract syntax tree more efficiently, general clean-up -- alphabetization
! etc.
!
! Revision 2.6  2000/11/30 00:23:58  pwagner
! functions properly moved here from Fill
!
! Revision 2.5  2000/11/29 17:35:30  pwagner
! Compiles now
!
! Revision 2.4  2000/11/29 00:27:54  pwagner
! Began changes to open old l2gp
!
! Revision 2.3  2000/11/16 01:02:16  vsnyder
! Correct an error message.
!
! Revision 2.2  2000/09/11 19:48:01  ahanzel
! Removed old log entries in file.
!
! Revision 2.1  2000/09/08 22:55:56  vsnyder
! Revised to use the tree output by the parser
!

