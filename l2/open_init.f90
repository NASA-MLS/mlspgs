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
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: FileNameLen, L1BInfo_T, TAI93_Range_T
!  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
!    &                         MLSMSG_Error, MLSMSG_FileOpen!, MLSMSG_Info
  use MLSPCF2, only: MLSPCF_L1B_OA_START, MLSPCF_L1B_RAD_END, &
    &                MLSPCF_L1B_RAD_START, MLSPCF_NOMEN_START, &
    &                mlspcf_pcf_start
  use MoreTree, only: Get_Spec_ID
  USE output_m, only: output
  USE PCFHdr, only: CreatePCFAnnotation
  use SDPToolkit, only: PGS_IO_Gen_closeF, PGS_IO_Gen_openF, &
    &                   Pgs_pc_getReference, PGS_S_SUCCESS, &
    &                   PGSd_IO_Gen_RSeqFrm, PGSTD_E_NO_LEAP_SECS
  use String_Table, only: Get_String !, L2CFUnit => INUNIT
  use Toggles, only: Gen, Levels, Switches, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, &
    &             SUB_ROSA, SUBTREE, DUMP_TREE_NODE, SOURCE_REF
  use TREE_TYPES, only: N_NAMED!, N_DOT
  use WriteMetadata, only: PCFData_T

  implicit none
  private
  public :: DestroyL1BInfo, OpenAndInitialize

  ! -----     Private declarations     ---------------------------------

  private :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  character(len=130) :: id = &
    "$id: open_init.f90,v 1.11 2000/06/19 22:40:51 lungu Exp $"
  character(len=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

    integer, parameter :: CCSDSLen=27
  integer, private :: ERROR

contains ! =====     Public Procedures     =============================

  ! ---------------------------------------------  DestroyL1BInfo  -----
  subroutine DestroyL1BInfo ( L1BInfo )
    type (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset
    integer :: STATUS ! from deallocate
    error = 0
    if (associated(l1bInfo%L1BRADIDs)) then
       deallocate( l1bInfo%L1BRADIDs, stat=status )
       if ( status /= 0 ) then
		 !	call MLSMessage ( MLSMSG_Error, ModuleName, &
       !     & MLSMSG_DeAllocate // "l1bInfo" )
		call announce_error(0, 'Error deallocating L1BRADIDs')
    endif
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
	 integer, external :: Pgs_pc_getFileSize

    !Local Variables
    logical, parameter :: DEBUG = .FALSE.
    integer, parameter :: CCSDSEndId = 10412
    integer, parameter :: CCSDSStartId = 10411

    character(len=*), parameter :: DEFAULTANTEXT= &
	 & 'PCF file number missing from PCF--add this line'
    character(len=CCSDSlen) CCSDSEndTime
    character(len=CCSDSlen) CCSDSStartTime
    integer :: ifl1
    integer :: L1FileHandle, L1_Version
    character (LEN=FileNameLen) :: L1physicalFilename
    integer :: returnStatus
    integer :: sd_id
    integer :: STATUS, size ! From allocate

      CHARACTER (LEN=FileNameLen) :: name

      INTEGER ::  indx, mlspcf_log, version

    integer :: pgs_td_utctotai, pgs_pc_getconfigdata

    error = 0
    if ( toggle(gen) ) call trace_begin ( "OpenAndInitialize" )

	if(DEBUG) call announce_error(0, &
	& 'Read the PCF into an annotation for file headers')

! Read the PCF into an annotation for file headers

      version = 1
      Status = Pgs_pc_getFileSize(mlspcf_pcf_start, version, size)
		if(Status == PGS_S_SUCCESS) then
	      CALL CreatePCFAnnotation(mlspcf_pcf_start, anText)
		else
			call announce_error(0, &
			& DEFAULTANTEXT)
			size = LEN(DEFAULTANTEXT) + 1
     		 ALLOCATE(anText(size), STAT=Status)
			anText(1:size-1) = DEFAULTANTEXT(1:size-1)
		endif

    ifl1 = 0

	if(DEBUG) call announce_error(0, 'Opening L1 RAD files')

    ! Open L1 RAD files
    do L1FileHandle = mlspcf_l1b_rad_start, mlspcf_l1b_rad_end

    ! Get the l1 file name from the PCF
    L1_Version = 1

      returnStatus = Pgs_pc_getReference(L1FileHandle, L1_Version, &
        & L1physicalFilename)

      if (returnStatus == PGS_S_SUCCESS) then

        ! Open the HDF file and initialize the SD interface

        ! Allocate L1BRADIDs

        allocate( l1bInfo%L1BRADIDs(10), stat=status )
        if ( status /= 0 ) then
!		  	call MLSMessage ( MLSMSG_Error, ModuleName, &
!          & MLSMSG_Allocate // "l1bInfo" )
			call announce_error(0, 'Allocation failed for L1BRADIDs')
			endif

        sd_id = sfstart(L1physicalFilename, DFACC_READ)
        if (sd_id == -1) then
!          call MLSMessage ( MLSMSG_Error, ModuleName, &
!            & "Error opening L1RAD file "//L1physicalFilename )
				call announce_error(0, &
				& 'Error opening L1RAD file: ' //L1physicalFilename)
        else
          ifl1 = ifl1 + 1
          l1bInfo%L1BRADIDs(ifl1) = sd_id
        end if
      end if
    end do ! L1FileHandle = mlspcf_l1b_rad_start, mlspcf_l1b_rad_end

    !if (ifl1 == 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
    !  & "Could not find any L1BRAD files" )

	if(DEBUG) call announce_error(0, 'Opening LOA file')
    ! Open L1OA File

    L1_Version = 1
    returnStatus = Pgs_pc_getReference(mlspcf_l1b_oa_start, L1_Version, &
      & L1physicalFilename)

    if (returnStatus == PGS_S_SUCCESS) then

      ! Open the HDF file and initialize the SD interface

      sd_id = sfstart(L1physicalFilename, DFACC_READ)
      if (sd_id == -1) then

!        call MLSMessage ( MLSMSG_Error, ModuleName, &
!          & "Error opening L1OA file "//L1physicalFilename )
			call announce_error(0, "Error opening L1OA file "//L1physicalFilename )
      else
        l1bInfo%L1BOAID = sd_id
      end if

    else

!      call MLSMessage ( MLSMSG_Error, ModuleName, "Could not find L1BOA file" )
			call announce_error(0, "Could not find L1BOA file" )

    end if

    ! Get the Start and End Times from PCF

    returnStatus = pgs_pc_getconfigdata (CCSDSStartId, CCSDSStartTime)
    if ( returnstatus /= PGS_S_SUCCESS ) then
	! call MLSMessage ( MLSMSG_Error, &
   !   & ModuleName, "Could not get CCSDS Start Time" )
			call announce_error(0, "Could not get CCSDS Start Time" )
	endif

    returnStatus = pgs_td_utctotai (CCSDSStartTime, processingrange%starttime)
    !   ??? Is PGSTD_E_NO_LEAP_SECS an OK status ???
    if ( returnstatus /= PGS_S_SUCCESS .and. &
      &  returnstatus /= PGSTD_E_NO_LEAP_SECS ) then
	!	call MLSMessage ( MLSMSG_Error, &
   !   & ModuleName, "Could not convert UTC Start time to TAI" )
			call announce_error(0, "Could not convert UTC Start time to TAI" )
	endif

    returnStatus = pgs_pc_getconfigdata (CCSDSEndId, CCSDSEndTime)
    if ( returnstatus /= PGS_S_SUCCESS ) then
	 !	call MLSMessage ( MLSMSG_Error, &
    !  & ModuleName, "Could not get CCSDS End Time" )
			call announce_error(0, "Could not get CCSDS End Time" )
	endif

    returnStatus = pgs_td_utctotai (CCSDSEndTime, processingrange%endtime)
    !   ??? Is PGSTD_E_NO_LEAP_SECS an OK status ???
    if ( returnstatus /= PGS_S_SUCCESS .and. &
      & returnstatus /= PGSTD_E_NO_LEAP_SECS) then
	!	call MLSMessage ( MLSMSG_Error, &
   !   & ModuleName, "Could not convert UTC Start time to TAI" )
			call announce_error(0, "Could not convert UTC End time to TAI" )
	endif

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
      IF (returnStatus /= PGS_S_SUCCESS) then
		!	CALL MLSMessage(MLSMSG_Error, &
      !                  ModuleName, 'Error retrieving log file name from PCF.')
			call announce_error(0, "Error retrieving log file name from PCF" )
	endif

      indx = INDEX(name, '/', .TRUE.)
      l2pcf%logGranID = name(indx+1:)
 
    if ( toggle(gen) ) then
      if ( levels(gen) > 0 .or. index(switches,'C') /= 0 ) &
        & call dump_L1B_database(ifl1, l1binfo, l2pcf, &
  			& CCSDSEndTime, CCSDSStartTime)
      call trace_end ( "OpenAndInit" )
    end if
    return

  end subroutine OpenAndInitialize

  ! -------------------------------------------------  dump_L1B_database  -----
  subroutine dump_L1B_database(num_l1b_files, l1binfo, l2pcf, &
  & CCSDSEndTime, CCSDSStartTime)
  
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

  call output('L1B database:', advance='yes')
  
  do i=1, num_l1b_files
    returnStatus = Pgs_pc_getReference(l1bInfo%L1BRADIDs(i), version, &
        & physicalFilename)
  	call output('fileid:   ')
	call output(l1bInfo%L1BRADIDs(i), advance='yes')
  	call output('name:   ')
  	call output(TRIM(physicalFilename), advance='yes')
  enddo

  call output('L1OA file:', advance='yes')
  
    returnStatus = Pgs_pc_getReference(l1bInfo%L1BOAID, version, &
        & physicalFilename)
  	call output('fileid:   ')
	call output(l1bInfo%L1BOAID, advance='yes')
  	call output('name:   ')
  	call output(TRIM(physicalFilename), advance='yes')

  	call output('Start Time:   ')
	call output(CCSDSStartTime, advance='yes')

  	call output('End Time:   ')
	call output(CCSDSEndTime, advance='yes')

  	call output('Output version:   ')
	call output(l2pcf%outputVersion, advance='yes')

  	call output('cycle:   ')
	call output(l2pcf%cycle, advance='yes')

  	call output('Log file name:   ')
  	call output(TRIM(l2pcf%logGranID), advance='yes')

  end subroutine dump_L1B_database

  ! ------------------------------------------------  announce_error  -----
  subroutine announce_error ( lcf_where, full_message, use_toolkit, &
  & error_number )
  
   ! Arguments
	
	integer, intent(in)    :: lcf_where
	character(LEN=*), intent(in)    :: full_message
	logical, intent(in), optional :: use_toolkit
	integer, intent(in), optional    :: error_number

	! Local
  logical :: just_print_it
  logical, parameter :: default_output_by_toolkit = .true.
	
	if(present(use_toolkit)) then
		just_print_it = use_toolkit
	elseif(default_output_by_toolkit) then
		just_print_it = .false.
	else
		just_print_it = .true.
	endif
	
	if(.not. just_print_it) then
    error = max(error,1)
    call output ( '***** At ' )

	if(lcf_where > 0) then
	    call print_source ( source_ref(lcf_where) )
		else
    call output ( '(no lcf node available)' )
		endif

    call output ( ': ' )
    call output ( "The " );
	if(lcf_where > 0) then
    call dump_tree_node ( lcf_where, 0 )
		else
    call output ( '(no lcf tree available)' )
		endif

		CALL output("Caused the following error:", advance='yes', &
		& from_where=ModuleName)
		CALL output(trim(full_message), advance='yes', &
		& from_where=ModuleName)
		if(present(error_number)) then
			CALL output('error number ', advance='no')
			CALL output(error_number, places=9, advance='yes')
		endif
	else
		print*, '***Error in module ', ModuleName
		print*, trim(full_message)
		if(present(error_number)) then
			print*, 'error number ', error_number
		endif
	endif

!===========================
  end subroutine announce_error
!===========================

end module Open_Init

!=============================================================================

!
! $Log$
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

