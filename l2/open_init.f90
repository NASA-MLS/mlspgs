! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Open_Init
!=============================================================================

! Opens and closes several files
! Creates and destroys the L1BInfo database

  use Hdf, only: DFACC_READ, SFSTART
  use Hdfeos, only: swopen, swclose
  use L2AUXData, only: L2AUXData_T, AddL2AUXToDatabase	!, ReadL2AUXData
  use L2GPData, only: L2GPData_T, AddL2GPToDatabase, ReadL2GPData
  use MLSCommon, only: L1BInfo_T, TAI93_Range_T
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error, MLSMSG_FileOpen
  use MLSPCF, only: MLSPCF_L1B_OA_START, MLSPCF_L1B_RAD_END, &
    &               MLSPCF_L1B_RAD_START, MLSPCF_NOMEN_START
  use MLSSignalNomenclature, only: ReadSignalsDatabase
  use MLSStrings, only: lowercase
  use SDPToolkit, only: PGS_IO_Gen_closeF, PGS_IO_Gen_openF, &
    &                   Pgs_pc_getReference, PGS_S_SUCCESS, &
    &                   PGSd_IO_Gen_RSeqFrm, PGSTD_E_NO_LEAP_SECS
  use String_Table, only: L2CFUnit => INUNIT
  use TREE, only: DECORATE, DECORATION, DUMP_TREE_NODE, NODE_ID, NSONS, &
    & SOURCE_REF, SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_NAMED, N_DOT
  use VectorsModule, only: AddVectorToDatabase, CreateVector, Dump, Vector_T, &
    & VectorTemplate_T

  implicit none
  private
  public :: CloseMLSCF, DestroyL1BInfo, OpenAndInitialize, OpenMLSCF, read_apriori

  ! -----     Private declarations     ---------------------------------

  private :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  character(len=130) :: id = &
     "$id: open_init.f90,v 1.11 2000/06/19 22:40:51 lungu Exp $"
  character(len=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  integer, parameter :: s_l2gp = 0   ! to be replaced by entry in init_tables_module
  integer, parameter :: s_l2aux = s_l2gp+1   ! to be replaced by entry in init_tables_module

contains ! =====     Public Procedures     =============================

  ! -------------------------------------------------  CloseMLSCF  -----
  subroutine CloseMLSCF

    character (len=32) :: Mnemonic
    character (len=256) :: Msg
    integer :: ReturnStatus

!   returnStatus = Pgs_io_gen_closeF ( L2CFUnit )

!   if ( returnStatus /= PGS_S_SUCCESS ) then
!     call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
!     call MLSMessage ( MLSMSG_Error, ModuleName, &
!       & 'Error closing L2CF:  '//mnemonic//' '//msg)
!   end if

  end subroutine CloseMLSCF

  ! ---------------------------------------------  DestroyL1BInfo  -----
  subroutine DestroyL1BInfo ( L1BInfo )
    type (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset
    integer :: STATUS ! from deallocate
    deallocate( l1bInfo%L1BRADIDs, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate // "l1bInfo" )
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
            & "Error opening L1RAD file"//L1physicalFilename )

        else
          ifl1 = ifl1 + 1
          l1bInfo%L1BRADIDs(ifl1) = sd_id
        end if
      end if
    end do ! L1FileHandle = mlspcf_l1b_rad_start, mlspcf_l1b_rad_end

    if (ifl1 == 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Could not find any L1BRAD files" )

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

    !processingrange%starttime = 220838406.65045166
    !processingrange%endtime = 220841983.31690979
    !processingrange%endtime = 220924616.64530945

    return

  end subroutine OpenAndInitialize

  ! --------------------------------------------------  OpenMLSCF  -----
  subroutine OpenMLSCF

    integer :: L2CF_Version
    character (len=32) :: Mnemonic
    character (len=256) :: Msg
    integer :: ReturnStatus

!   Open the MLSCF as a generic file for reading
    L2CF_Version = 1
!   returnStatus = Pgs_io_gen_openF ( mlspcf_l2cf_start, PGSd_IO_Gen_RSeqFrm, &
!     &                               0, L2CFUnit, L2CF_Version)

!   if ( returnStatus /= PGS_S_SUCCESS ) then

!     call Pgs_smf_getMsg ( returnStatus, mnemonic, msg )
!     call MLSMessage ( MLSMSG_Error, ModuleName, &
!       & "Error opening MLSCF:  "//mnemonic//"  "//msg)
!   end if

  end subroutine OpenMLSCF

  ! --------------------------------------------------  read_apriori  -----
  subroutine read_apriori ( root, L2GPDatabase, l2auxDatabase)

	! Read in  a priori data from l2gp files and l2aux files

    ! Dummy arguments
    integer, intent(in) :: ROOT    ! Of the Read a priori section in the AST
    TYPE (L2GPData_T), DIMENSION(:), POINTER :: L2GPDatabase
    type (L2AUXData_T), dimension(:), pointer :: l2auxDatabase

    !Local Variables
    integer :: I, J                ! Loop indices for section, spec
    integer :: KEY                 ! Definitely n_named
    type (Vector_T) :: newVector
    integer :: SON                 ! Of root, an n_spec_args or a n_named
    integer :: templateIndex       ! In the template database
    integer :: vectorIndex         ! In the vector database
    integer :: vectorName          ! Sub-rosa index
    integer :: quantityName        ! Sub-rosa index
    integer :: sourceName          ! Sub-rosa index
    integer :: FileType            ! Sub-rosa index of either 'l2gp' or 'l2aux'
    integer :: fileHandle          ! fileHandle of a priori data file
    integer :: NumProfs            ! number of profiles actually read
    CHARACTER (LEN=24) :: FileTypeString     ! 'l2gp' or 'l2aux'
    CHARACTER (LEN=480) :: msr     ! Error message if can't find file
    integer :: JName               ! index to spec name of 3rd son of root
    CHARACTER*16 :: JString        ! index to spec name of 3rd son of root
    integer :: FileRoot            ! index to spec name of file field
    integer :: SourceRoot          ! index to spec name of source field
    integer :: FileName            ! subrosa indexed of name in file='name'
    CHARACTER*16 :: FileNameString ! actual literal file name
    CHARACTER*16 :: SourceNameString ! actual literal source name
    TYPE (L2GPData_T) :: L2GP
    TYPE (L2AUXData_T) :: L2AUX

! Toolkit Functions

      INTEGER, EXTERNAL :: swclose, swopen

	! Assume specifications take the following form:
        !   vName: fileType, file='fileName', source='fieldName'
        ! Currently, fileType is restricted to one of
        !   'l2gp' or 'l2aux'
    do i = 2, nsons(root)-1 ! Skip the section name at begin and end
      son = subtree(i,root)
      if ( node_id(son) == n_named ) then ! Is spec labeled?
        key = subtree(2,son)
        vectorName = sub_rosa(subtree(1,son))
      else
        key = son
        vectorName = 0
      end if

      ! Node_id(key) is now n_spec_args.

        if ( nsons(key) /= 4 ) call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Wrong number of fields in read a priori')

        FileType = sub_rosa(subtree(2,son))
        CALL get_string(FileType, FileTypeString)
	FileTypeString = lowercase(FileTypeString)

        ! Now parse file and field names
        ! Does file spec precede field spec?
        J = subtree(3, key)
        JName = decoration(subtree(1,decoration(subtree(1,J))))
        CALL get_string(JName, JString)
        JString = lowercase(JString)
        IF(JString.EQ.'file') THEN
 		! Yes
        	fileRoot = subtree(3, key)
        	sourceRoot = subtree(4, key)
	ELSE
 		! No
        	fileRoot = subtree(4, key)
        	sourceRoot = subtree(3, key)
	ENDIF
	FileName=sub_rosa(subtree(2, fileRoot))
        CALL get_string(FileName, FileNameString)
	sourceName=sub_rosa(subtree(2, sourceRoot))
        CALL get_string(sourceName, sourceNameString)

            fileHandle = swopen(FileNameString, DFACC_READ)
            IF (fileHandle == -1) THEN
               msr = MLSMSG_Fileopen // FileNameString
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

      select case( FileTypeString )
      case ( 'l2gp' )

        vectorIndex = decoration(decoration(subtree(2,subtree(1,key))))

        ! Create the l2gp, and add it to the database.
        CALL SetupNewL2GPRecord ( l2gp )

        call decorate ( key, AddL2GPToDatabase( L2GPDatabase, l2gp ) )

        ! That's the end of the create operation
        

           CALL ReadL2GPData(fileHandle, sourceNameString, L2GPDatabase(vectorIndex),&
           & numProfs)


      case ( 'l2aux' )

        vectorIndex = decoration(decoration(subtree(2,subtree(1,key))))

        ! Create the l2aux, and add it to the database.
        CALL SetupNewL2AUXRecord ( l2aux )

        call decorate ( key, AddL2AUXToDatabase( L2AUXDatabase, l2aux ) )

        ! That's the end of the create operation
        
! Need to add this routine to L2AUXData.f90 before uncommenting this line
!           CALL ReadL2AUXData(fileHandle, sourceNameString, L2GPDatabase(vectorIndex),&
!           & numProfs)


      case default ! Can't get here if tree_checker worked correctly
      end select

            fileHandle = swclose(fileHandle)
            IF (fileHandle == -1) THEN
               msr = 'Failed to close file ' // FileNameString
               CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
            ENDIF

    end do

  end subroutine read_apriori

!=============================================================================
end module Open_Init
!=============================================================================

!
! $Log$
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

