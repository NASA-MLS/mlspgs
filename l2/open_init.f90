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
    &                MLSPCF_L1B_RAD_START, MLSPCF_NOMEN_START, &
    &                mlspcf_pcf_start
  use MoreTree, only: Get_Spec_ID
  USE PCFHdr, only: CreatePCFAnnotation
  use SDPToolkit, only: PGS_IO_Gen_closeF, PGS_IO_Gen_openF, &
    &                   Pgs_pc_getReference, PGS_S_SUCCESS, &
    &                   PGSd_IO_Gen_RSeqFrm, PGSTD_E_NO_LEAP_SECS
  use String_Table, only: Get_String !, L2CFUnit => INUNIT
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, &
    &             SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_NAMED!, N_DOT
  use WriteMetadata, only: PCFData_T

  implicit none
  private
  public :: Close_MLSCF, DestroyL1BInfo, OpenAndInitialize, Open_MLSCF

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
  subroutine OpenAndInitialize ( processingRange, l1bInfo, l2pcf, anText )

    ! Opens L1 RAD files
    ! Opens L1OA file
    ! Opens and reads the signal nomenclature file
    ! Gets the start and end times from the PCF

    ! Arguments

    type (TAI93_Range_T) :: processingRange ! Data processing range
    type (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset
	 type(PCFData_T) :: l2pcf
    CHARACTER (LEN=1), POINTER :: anText(:)

    !Local Variables
    integer, parameter :: CCSDSEndId = 10412
    integer, parameter :: CCSDSLen=27
    integer, parameter :: CCSDSStartId = 10411

    character(len=CCSDSlen) CCSDSEndTime
    character(len=CCSDSlen) CCSDSStartTime
    integer :: ifl1
    integer :: L1FileHandle, L1_Version
    character (LEN=132) :: L1physicalFilename
    integer :: returnStatus
    integer :: sd_id
    integer :: STATUS ! From allocate

      CHARACTER (LEN=FileNameLen) :: name

      INTEGER ::  indx, mlspcf_log, version

    integer :: pgs_td_utctotai, pgs_pc_getconfigdata

! Read the PCF into an annotation for file headers

      CALL CreatePCFAnnotation(mlspcf_pcf_start, anText)

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

	! -- May need to change these if we have
	!    configuration parameters introduced into PCF file
	!    for now -- outputversion and cycle are hard-wired below
	! --

	! Here's where we define the components of l2pcf

!      returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_OutputVersion, &
!                                          l2pcf%outputVersion)
	l2pcf%outputVersion = '1'
	
!      returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_Cycle, l2pcf%cycle)
	l2pcf%cycle = '1'
	
	l2pcf%startutc = CCSDSStartTime
	l2pcf%endutc = CCSDSEndTime

! Get the name of the log file from the PCF

      version = 1
      mlspcf_log = 10101			! This seems to be hard-wired into PCF

      returnStatus = Pgs_pc_getReference(mlspcf_log, version, name)
      IF (returnStatus /= PGS_S_SUCCESS) CALL MLSMessage(MLSMSG_Error, &
                        ModuleName, 'Error retrieving log file name from PCF.')
      indx = INDEX(name, '/', .TRUE.)
      l2pcf%logGranID = name(indx+1:)
 
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

end module Open_Init
!=============================================================================

!
! $Log$
! Revision 2.26  2001/04/03 20:51:27  pwagner
! Added anText; deleted read_apriori
!
! Revision 2.25  2001/04/02 23:39:09  pwagner
! Now fills components of l2pcf
!
! Revision 2.24  2001/03/28 19:07:59  vsnyder
! Finish removing use of MLSSignalNomenclature
!
! Revision 2.23  2001/03/28 19:07:08  vsnyder
! Remove use of MLSSignalNomenclature
!
! Revision 2.22  2001/03/15 21:18:57  vsnyder
! Use Get_Spec_ID instead of decoration(subtree...
!
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

