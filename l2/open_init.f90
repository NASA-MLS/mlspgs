! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module Open_Init

  ! Opens and closes several files
  ! Creates and destroys the L1BInfo database

  use Hdf, only: DFACC_READ, SFSTART, SFEND
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: FileNameLen, L1BInfo_T, TAI93_Range_T
  use MLSL2Options, only: PUNISH_FOR_INVALID_PCF, PUNISH_FOR_NO_L1BRAD, &
  &                        PUNISH_FOR_NO_L1BOA, PENALTY_FOR_NO_METADATA
  use MLSMessageModule, only: MLSMessage, &
    &                         MLSMSG_Error!, MLSMSG_FileOpen, MLSMSG_Info
  use MLSPCF2, only: MLSPCF_L1B_OA_START, MLSPCF_L1B_RAD_END, &
    &                MLSPCF_L1B_RAD_START, &
    &                mlspcf_l2_param_InputVersion, &
    &                mlspcf_l2_param_PGEVersion, &
    &                mlspcf_l2_param_Cycle, &
    &                mlspcf_l2_param_CCSDSStartId, &
    &                mlspcf_l2_param_CCSDSEndId, &
    &                mlspcf_l2_param_spec_keys, &
    &                mlspcf_l2_param_spec_hash, &
    &                mlspcf_pcf_start
  use MLSStrings, only: LowerCase
  use Output_m, only: Output
  use PCFHdr, only: CreatePCFAnnotation
  use SDPToolkit, only: &
    &                   Pgs_pc_getReference, PGS_S_SUCCESS, &
    &                   PGSTD_E_NO_LEAP_SECS
  use Toggles, only: Gen, Levels, Switches, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use TREE, only: DUMP_TREE_NODE, SOURCE_REF
  use WriteMetadata, only: PCFData_T, MCFCASESENSITIVE

  implicit none
  private
  public :: DestroyL1BInfo, OpenAndInitialize

  ! -----     Private declarations     ---------------------------------

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = &
    "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  integer, parameter :: CCSDSLen=27
  integer, parameter :: ILLEGALL1BRADID=-1      ! something sfstart should catch
  integer, parameter :: MAXNUML1BRADIDS=10      ! In case more than one
  integer, private :: ERROR
  
contains ! =====     Public Procedures     =============================

  ! ---------------------------------------------  DestroyL1BInfo  -----
  subroutine DestroyL1BInfo ( L1BInfo )
    type (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset
    integer :: STATUS ! from deallocate
    integer :: id
    error = 0
    if ( toggle(gen) ) call trace_begin ( "DESTROYL1BInfo" )

    if ( associated(l1bInfo%L1BRADIDs) ) then
      do id=1, SIZE(l1bInfo%L1BRADIDs)
         if(l1bInfo%L1BRADIDs(id) /= ILLEGALL1BRADID) then
          STATUS = sfend(l1bInfo%L1BRADIDs(id))
         endif
      enddo
      deallocate( l1bInfo%L1BRADIDs, stat=status )
      if ( status /= 0 ) then
        ! call MLSMessage ( MLSMSG_Error, ModuleName, &
        !   & MLSMSG_DeAllocate // "l1bInfo" )
        call announce_error ( 0, 'Error deallocating L1BRADIDs' )
      end if
    end if
    
    if(l1bInfo%L1BOAID /= ILLEGALL1BRADID) then
      STATUS = sfend(l1bInfo%L1BOAID)
    endif

    if ( toggle(gen) ) then
      call trace_end ( "DESTROYL1BInfo" )
    end if
  end subroutine DestroyL1BInfo


  ! ------------------------------------------  OpenAndInitialize  -----
  subroutine OpenAndInitialize ( processingRange, l1bInfo, l2pcf )

    ! Opens L1 RAD files
    ! Opens L1OA file
    ! Opens and reads the signal nomenclature file
    ! Gets the start and end times from the PCF

    ! Arguments

    type (TAI93_Range_T) :: processingRange ! Data processing range
    type (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset
    type(PCFData_T) :: l2pcf
    integer, external :: Pgs_pc_getFileSize

    !Local Variables
    logical, parameter :: DEBUG = .FALSE.
    integer, parameter :: CCSDSEndId = 10412    ! Illegal PCFid
    integer, parameter :: CCSDSStartId = 10411    ! Illegal PCFid

   ! The following parameters will only be needed if PCF ids are missing
    character(len=*), parameter :: DEFAULTANTEXT= &
      & 'PCF file number missing from PCF--add this line'
    character(len=*), parameter :: DEFAULT_SPEC_KEYS= &
      & 'temp,gph,h2o,hno3,o3,hcl,clo,co,n2o,oh,rhi,so2,ho2,bro,hocl,hcn,cirrus-ice,others'
    character(len=*), parameter :: DEFAULT_SPEC_HASH= &
      & 't,z,h2o,hno3,o3,hcl,clo,co,n2o,oh,rhi,so2,ho2,bro,hocl,hcn,ice,oth'

    character(len=CCSDSlen) CCSDSEndTime
    character(len=CCSDSlen) CCSDSStartTime
    integer :: Ifl1
    integer :: L1FileHandle, L1_Version
    character (len=fileNameLen) :: L1physicalFilename
    integer :: ReturnStatus
    integer :: Sd_id
    integer :: Status, Size ! From allocate

    character (len=fileNameLen) :: name

    integer :: Indx, Mlspcf_log, Version

    integer :: pgs_td_utctotai, pgs_pc_getconfigdata

    error = 0
    if ( toggle(gen) ) call trace_begin ( "OpenAndInitialize" )

! Read the PCF into an annotation for file headers

    version = 1
    Status = Pgs_pc_getFileSize(mlspcf_pcf_start, version, size)
    if ( Status == PGS_S_SUCCESS ) then
      call createPCFAnnotation(mlspcf_pcf_start, l2pcf%anText)
    else
      call announce_error ( 0, DEFAULTANTEXT )
      size = LEN(DEFAULTANTEXT) + 1
      allocate ( l2pcf%anText(size), STAT=Status )
      l2pcf%anText(1:size-1) = DEFAULTANTEXT(1:size-1)
      error = PENALTY_FOR_NO_METADATA
    end if

    ifl1 = 0

    ! Open L1 RAD files
    do L1FileHandle = mlspcf_l1b_rad_start, mlspcf_l1b_rad_end

    ! Get the l1 file name from the PCF
    L1_Version = 1

      returnStatus = Pgs_pc_getReference(L1FileHandle, L1_Version, &
        & L1physicalFilename)

      if ( returnStatus == PGS_S_SUCCESS ) then

        ! Open the HDF file and initialize the SD interface

        ! Allocate L1BRADIDs, initialize them to illegal values

      if(.NOT. associated(l1bInfo%L1BRADIDs)) then
        allocate ( l1bInfo%L1BRADIDs(MAXNUML1BRADIDS), stat=status )
        l1bInfo%L1BRADIDs = ILLEGALL1BRADID
        if ( status /= 0 ) &
          & call announce_error ( 0, 'Allocation failed for L1BRADIDs' )
      endif
        sd_id = sfstart(L1physicalFilename, DFACC_READ)
        if ( sd_id == -1 ) then
          call announce_error ( 0, &
            & 'Error opening L1RAD file: ' //L1physicalFilename)
        else
          ifl1 = ifl1 + 1
          l1bInfo%L1BRADIDs(ifl1) = sd_id
        end if
      end if
    end do ! L1FileHandle = mlspcf_l1b_rad_start, mlspcf_l1b_rad_end

    if ( ifl1 == 0 .AND. PUNISH_FOR_NO_L1BRAD ) &
      &  call announce_error ( 0, "Could not find any L1BRAD files" )

    ! Open L1OA File

    l1bInfo%L1BOAID = ILLEGALL1BRADID
    L1_Version = 1
    returnStatus = Pgs_pc_getReference(mlspcf_l1b_oa_start, L1_Version, &
      & L1physicalFilename)

    if ( returnStatus == PGS_S_SUCCESS ) then

      ! Open the HDF file and initialize the SD interface

      sd_id = sfstart(L1physicalFilename, DFACC_READ)
      if ( sd_id == -1 ) then

!        call MLSMessage ( MLSMSG_Error, ModuleName, &
!          & "Error opening L1OA file "//L1physicalFilename )
        call announce_error ( 0, "Error opening L1OA file "//L1physicalFilename )
      else
        l1bInfo%L1BOAID = sd_id
      end if

    else if ( PUNISH_FOR_NO_L1BOA ) then
      call announce_error ( 0, "Could not find L1BOA file" )
    end if

    ! Get the Start and End Times from PCF
    ! Temporarily we allow the use of older PCFids: CCSDSStartId, CCSDEndId

    returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_CCSDSStartId, &
                                           CCSDSStartTime)
    if ( returnstatus /= PGS_S_SUCCESS .and. PUNISH_FOR_INVALID_PCF ) then
      call announce_error ( 0, "Missing pcf param: CCSDSStartTime" )
    else if ( returnstatus /= PGS_S_SUCCESS ) then
      returnStatus = pgs_pc_getconfigdata (CCSDSStartId, CCSDSStartTime)
      if ( returnstatus /= PGS_S_SUCCESS ) &
        & call announce_error ( 0, "Could not get CCSDS Start Time" )
    end if

    returnStatus = pgs_td_utctotai (CCSDSStartTime, processingrange%starttime)
    !   ??? Is PGSTD_E_NO_LEAP_SECS an OK status ???
    if ( returnstatus /= PGS_S_SUCCESS .and. &
      &  returnstatus /= PGSTD_E_NO_LEAP_SECS ) &
        & call announce_error ( 0, "Could not convert UTC Start time to TAI" )

   returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_CCSDSEndId, &
                                          CCSDSEndTime)
    if ( returnstatus /= PGS_S_SUCCESS .and. PUNISH_FOR_INVALID_PCF ) then
      call announce_error ( 0, "Missing pcf param: CCSDSEndTime" )
    else if ( returnstatus /= PGS_S_SUCCESS ) then
      returnStatus = pgs_pc_getconfigdata (CCSDSEndId, CCSDSEndTime)
      if ( returnstatus /= PGS_S_SUCCESS ) &
        & call announce_error ( 0, "Could not get CCSDS End Time" )
    end if

    returnStatus = pgs_td_utctotai (CCSDSEndTime, processingrange%endtime)
    !   ??? Is PGSTD_E_NO_LEAP_SECS an OK status ???
    if ( returnstatus /= PGS_S_SUCCESS .and. &
      & returnstatus /= PGSTD_E_NO_LEAP_SECS) &
        & call announce_error ( 0, "Could not convert UTC End time to TAI" )

    l2pcf%startutc = CCSDSStartTime
    l2pcf%endutc = CCSDSEndTime

    ! Here's where we define the non-time components of l2pcf

    returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_inputVersion, &
                                          l2pcf%inputVersion)
    if ( returnstatus /= PGS_S_SUCCESS .and. PUNISH_FOR_INVALID_PCF ) then
      call announce_error ( 0, "Missing pcf param: input version" )
    else if ( returnstatus /= PGS_S_SUCCESS ) then
      l2pcf%inputVersion = 'V0-5'
    end if

    returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_PGEVersion, &
                                          l2pcf%PGEVersion)
    if ( returnstatus /= PGS_S_SUCCESS .and. PUNISH_FOR_INVALID_PCF ) then
      call announce_error ( 0, "Missing pcf param: output version" )
    else if ( returnstatus /= PGS_S_SUCCESS ) then
      l2pcf%PGEVersion = 'V0-5'
    end if
	
    returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_Cycle, l2pcf%cycle)

    if ( returnstatus /= PGS_S_SUCCESS .and. PUNISH_FOR_INVALID_PCF ) then
      call announce_error ( 0, "Missing pcf param: cycle" )
    else if ( returnstatus /= PGS_S_SUCCESS ) then
      l2pcf%cycle = '1'
    end if
	
    returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_spec_keys, l2pcf%spec_keys)

    if ( returnstatus /= PGS_S_SUCCESS .and. PUNISH_FOR_INVALID_PCF ) then
      call announce_error ( 0, "Missing pcf param: spec_keys" )
    else if ( returnstatus /= PGS_S_SUCCESS ) then
      l2pcf%spec_keys = DEFAULT_SPEC_KEYS
    end if
	
    returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_spec_hash, l2pcf%spec_hash)

    if ( returnstatus /= PGS_S_SUCCESS .and. PUNISH_FOR_INVALID_PCF ) then
      call announce_error ( 0, "Missing pcf param: spec_hash" )
    else if ( returnstatus /= PGS_S_SUCCESS ) then
      l2pcf%spec_hash = DEFAULT_SPEC_HASH
    end if
	
    if ( .NOT. MCFCASESENSITIVE ) then
      l2pcf%spec_hash = LowerCase(l2pcf%spec_hash)
      l2pcf%spec_keys = LowerCase(l2pcf%spec_keys)
    end if

! Get the name of the log file from the PCF

    version = 1
    mlspcf_log = 10101			! This seems to be hard-wired into PCF

    returnStatus = Pgs_pc_getReference(mlspcf_log, version, name)
    if ( returnStatus /= PGS_S_SUCCESS ) &
      & call announce_error ( 0, "Error retrieving log file name from PCF" )

    indx = INDEX(name, '/', .TRUE.)
    l2pcf%logGranID = name(indx+1:)
 
    if ( error /= 0 ) &
      & call MLSMessage(MLSMSG_Error,ModuleName, &
        & 'Problem with open_init section')

    if ( toggle(gen) ) then
      if ( levels(gen) > 0 .or. index(switches,'O') /= 0 ) &
        & call dump_L1B_database ( ifl1, l1binfo, l2pcf, &
          & CCSDSEndTime, CCSDSStartTime )
      call trace_end ( "OpenAndInit" )
    end if
    return

  end subroutine OpenAndInitialize

  ! ------------------------------------------  Dump_L1B_database  -----
  subroutine Dump_L1B_database ( Num_l1b_files, L1binfo, L2pcf, &
    & CCSDSEndTime, CCSDSStartTime )
  
    ! Dump info obtained during OpenAndInitialize:
    ! L1B databse
    ! L1OA file
    ! Start and end times
    ! output version
    ! cycle number
    ! logfile name
  

    ! Arguments
    integer, intent(in) :: num_l1b_files
    type (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset
    type(PCFData_T) :: l2pcf
    character(len=CCSDSlen) CCSDSEndTime
    character(len=CCSDSlen) CCSDSStartTime
	
    ! Local

    character (LEN=FileNameLen) :: physicalFilename
    integer :: i, returnStatus, version

    ! Begin
    version = 1

    call output ( 'L1B database:', advance='yes' )
  
    if ( num_l1b_files > 0 ) then
      do i = 1, num_l1b_files
        returnStatus = Pgs_pc_getReference(l1bInfo%L1BRADIDs(i), version, &
        & physicalFilename)
  	call output ( 'fileid:   ' )
	call output ( l1bInfo%L1BRADIDs(i), advance='yes' )
  	call output ( 'name:   ' )
  	call output ( TRIM(physicalFilename), advance='yes' )
      end do
    else
      call output ( '(empty database)', advance='yes' )
    end if

    call output ( 'L1OA file:', advance='yes' )
  
    returnStatus = Pgs_pc_getReference(l1bInfo%L1BOAID, version, &
      & physicalFilename)
    if ( returnStatus == PGS_S_SUCCESS ) then
      call output ( 'fileid:   ' )
      call output ( l1bInfo%L1BOAID, advance='yes' )
      call output ( 'name:   ' )
      call output ( TRIM(physicalFilename), advance='yes' )
    else
      call output ( '(file unknown)', advance='yes' )
    end if

    call output ( 'Start Time:   ' )
    call output ( CCSDSStartTime, advance='yes' )

    call output ( 'End Time:   ' )
    call output ( CCSDSEndTime, advance='yes' )

    call output ( 'PGE version:   ' )
    call output ( l2pcf%PGEVersion, advance='yes' )

    call output ( 'input version:   ' )
    call output ( l2pcf%InputVersion, advance='yes' )

    call output ( 'cycle:   ' )
    call output ( l2pcf%cycle, advance='yes' )

    call output ( 'Log file name:   ' )
    call output ( TRIM(l2pcf%logGranID), advance='yes' )

    call output ( 'l2gp species name keys:   ' )
    call output ( TRIM(l2pcf%spec_keys), advance='yes' )

    call output ( 'corresponding mcf hash:   ' )
    call output ( TRIM(l2pcf%spec_hash), advance='yes' )

  end subroutine dump_L1B_database

  ! ---------------------------------------------  Announce_Error  -----
  subroutine Announce_Error ( Lcf_where, Full_message, Use_toolkit, &
    & Error_number )
  
    ! Arguments

    integer, intent(in) :: Lcf_where
    character(LEN=*), intent(in) :: Full_message
    logical, intent(in), optional :: Use_toolkit
    integer, intent(in), optional :: Error_number

    ! Local
    logical :: Just_print_it
    logical, parameter :: Default_output_by_toolkit = .true.

    just_print_it = .not. default_output_by_toolkit
    if ( present(use_toolkit) ) just_print_it = .not. use_toolkit

    if ( .not. just_print_it ) then
      error = max(error,1)
      call output ( '***** At ' )

      if ( lcf_where > 0 ) then
        call print_source ( source_ref(lcf_where) )
      else
        call output ( '(no lcf node available)' )
      end if

      call output ( ": The " );
      if ( lcf_where > 0 ) then
        call dump_tree_node ( lcf_where, 0 )
      else
        call output ( '(no lcf tree available)' )
      end if

      call output ( " Caused the following error:", advance='yes', &
       & from_where=ModuleName)
      call output ( trim(full_message), advance='yes', &
        & from_where=ModuleName)
      if ( present(error_number) ) then
        call output ( 'Error number ', advance='no' )
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
  end subroutine Announce_Error
!===========================

end module Open_Init

!=============================================================================

!
! $Log$
! Revision 2.42  2001/04/20 18:16:29  pwagner
! Added sfend on L1B files to DESTROYL1BInfo
!
! Revision 2.41  2001/04/19 23:51:40  pwagner
! Moved anText to become component of PCFData_T
!
! Revision 2.40  2001/04/18 17:33:29  vsnyder
! I forgot what I REALLY did, but I also made numerous cosmetic changes
!
! Revision 2.39  2001/04/16 23:45:59  pwagner
! Gets penalty settings from MLSL2Options
!
! Revision 2.38  2001/04/16 17:45:54  pwagner
! Makes l2pcf hash, keys lowercase if not MCFCASESENSITIVE
!
! Revision 2.37  2001/04/12 00:20:44  pwagner
! With new PCF params
!
! Revision 2.32  2001/04/10 23:00:29  pwagner
! Keeps track of whether to quit if no metadata
!
! Revision 2.31  2001/04/06 20:20:43  vsnyder
! Improve an error message
!
! Revision 2.30  2001/04/06 18:01:00  pwagner
! Checks on pcf number before PCFCreateAnnotation
!
! Revision 2.29  2001/04/05 23:44:53  pwagner
! Fixed tiny error
!
! Revision 2.28  2001/04/05 23:40:50  pwagner
! Deleted open_mlscf and close_mlscf and all MLSMessages
!
! Revision 2.27  2001/04/04 23:46:18  pwagner
! Added trace_*, dump_l1b_database
!
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

