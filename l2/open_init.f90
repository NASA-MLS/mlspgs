! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Open_Init

  ! Opens and closes several files
  ! Creates and destroys the L1BInfo database

  use HIGHOUTPUT, only: BEVERBOSE
  use MLSCOMMON, only: FILENAMELEN, MLSFILE_T, TAI93_RANGE_T
  use MLSL2OPTIONS, only: SPECIALDUMPFILE, TOOLKIT
  use MLSSTRINGLISTS, only: CATLISTS, GETSTRINGELEMENT, NUMSTRINGELEMENTS, &
    & SWITCHDETAIL
  use OUTPUT_M, only: BLANKS, OUTPUT

  implicit none

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (subroutines and functions)
! OpenAndInitialize               Gets run parameters from pcf
! === (end of toc) ===

  private
  public :: OpenAndInitialize

  ! -----     Private declarations     ---------------------------------

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  integer, parameter :: MLSPCF_LOG = 10101 ! This seems to be hard-wired into PCF
  integer, parameter :: CCSDSLen=27
  logical, parameter :: ALWAYSGETSWITCHESCONFIG = .false.
  integer, private ::   ERROR
  
contains ! =====     Public Procedures     =============================

  ! ------------------------------------------  OpenAndInitialize  -----
  subroutine OpenAndInitialize ( processingRange, filedatabase )

    ! Opens L1 RAD files
    ! Opens L1OA file
    ! Gets the start and end times from the PCF
    ! Gets other user parameters from PCF
    
    ! Note: toolkitless runs bypass all that, getting
    ! the necessary inputs via the l2cf and
    ! processing user inputs during global settings section

    use DATES_MODULE, only: UTC_TO_YYYYMMDD
    use HDF, only: DFACC_RDONLY
    use INTRINSIC, only: L_HDF
    use L1BDATA, only: FINDMAXMAF, READL1BATTRIBUTE
    use L2GPDATA, only: COL_SPECIES_KEYS, COL_SPECIES_HASH
    use MLSFILES, only: WILDCARDHDFVERSION, &
      & ADDFILETODATABASE, INITIALIZEMLSFILE, MLS_OPENFILE
    use MLSKINDS, only: R8
    use MLSL2TIMINGS, only: SECTION_TIMES, TOTAL_TIMES
    use MLSMESSAGEMODULE, only: MLSMSG_ERROR, MLSMSG_WARNING, MLSMESSAGE
    use MLSPCF2, only: MLSPCF_L1B_OA_START, MLSPCF_L1B_RAD_END, &
      &                MLSPCF_L1B_RAD_START, &
      &                MLSPCF_L2_PARAM_PGEVERSION, &
      &                MLSPCF_L2_PARAM_CYCLE, &
      &                MLSPCF_L2_PARAM_CCSDSSTARTID, &
      &                MLSPCF_L2_PARAM_CCSDSENDID, &
      &                MLSPCF_L2_PARAM_COL_SPEC_KEYS, &
      &                MLSPCF_L2_PARAM_COL_SPEC_HASH, &
      &                MLSPCF_L2_PARAM_SPEC_KEYS, &
      &                MLSPCF_L2_PARAM_SPEC_HASH, &
      &                MLSPCF_L2_PARAM_SWITCHES, &
      &                MLSPCF_PCF_START
    use MLSSTRINGS, only: LOWERCASE
    use PCFHDR, only: GLOBALATTRIBUTES, CREATEPCFANNOTATION, FILLTAI93ATTRIBUTE
    use SDPTOOLKIT, only: MAX_ORBITS, PGS_PC_GETFILESIZE, PGS_TD_UTCTOTAI,&
      &    PGS_PC_GETCONFIGDATA, PGS_PC_GETREFERENCE, PGS_S_SUCCESS, &
      &    PGSTD_E_NO_LEAP_SECS
    use TIME_M, only: TIME_NOW
    use TOGGLES, only: GEN, LEVELS, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use WRITEMETADATA, only: L2PCF, MCFCASESENSITIVE

    ! Arguments

    type (TAI93_Range_T) :: processingRange ! Data processing range
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE

    ! Local Variables
    ! logical, parameter :: DEBUG = .FALSE.
    logical, parameter :: CaseSensitive = .true.

    character(len=CCSDSlen)      :: CCSDSEndTime
    character(len=CCSDSlen)      :: CCSDSStartTime
    ! The following parameters will only be needed if PCF ids are missing
    character(len=*), parameter :: DEFAULTANTEXT= &
      & 'PCF file number missing from PCF--add this line'

    integer                      :: Details   ! How much info about l1b files to dump
    character(len=len(switches)) :: extra_switches

    integer                      :: Ifl1
    integer                      :: Indx, Version                            

    integer                      :: l1bFlag
    type (MLSFile_T)             :: L1BFile
    integer                      :: L1FileHandle
    integer :: Me = -1           ! String index for trace
    character (len=fileNameLen)  :: name
    integer                      :: numFiles

    integer                      :: ReturnStatus
    character(len=16)            :: RUNID
    integer                      :: Status, Size ! From allocate

    real                         :: T1, T2                      ! for timing 
    logical                      :: TIMING                                   
    ! Executable
    integer :: OrbNum(max_orbits) = 0
    real(r8) :: OrbPeriod(max_orbits) = 0.0

    ! Executable code
    call trace_begin ( me, "OpenAndInitialize", cond=toggle(gen) )
    timing = section_times
    if ( timing ) call time_now ( t1 )

    error = 0

! Read the PCF into an annotation for file headers

    if ( .not. TOOLKIT ) then
       Status = PGS_S_SUCCESS - 1
    else
     version = 1
     Status = Pgs_pc_getFileSize(mlspcf_pcf_start, version, size)
    end if
   
    if ( Status == PGS_S_SUCCESS ) then
      call createPCFAnnotation(mlspcf_pcf_start, l2pcf%anText)
    else
      if ( TOOLKIT ) then
        call announce_error ( DEFAULTANTEXT )
        error = 1
      endif
      size = LEN(DEFAULTANTEXT) + 1
      allocate ( l2pcf%anText(size), STAT=Status )
      l2pcf%anText(1:size-1) = DEFAULTANTEXT(1:size-1)
    end if

   ! Initialize run parameters: unless reset, either by pcf or l2cf,
   ! these determine how to run
   GlobalAttributes%ProcessLevel = '2'
   CCSDSStartTime = '(undefined)'
   CCSDSEndTime = '(undefined)'
   l2pcf%startutc = '(undefined)'
   l2pcf%endutc = '(undefined)'
   l2pcf%cycle = ' '
   l2pcf%PGEVersion = ' '
   l2pcf%logGranID = '(not applicable)'      ! will not create a Log file
   l2pcf%spec_keys = '(not applicable)'      ! will not create metadata
   l2pcf%spec_hash = '(not applicable)'      ! will not create metadata

   if( .not. TOOLKIT ) then
     if ( levels(gen) > 0 .or. switchDetail(switches,'pcf') > -1 ) then
       call output('====== No parameters or radiances read :: no pcf ======', &
         & advance='yes')
       call output('====  (These must be supplied through the l2cf) ======', &
         & advance='yes')
     end if
     call trace_end ( "OpenAndInitialize", cond=toggle(gen) )
     if ( timing ) call sayTime
     return
   end if

    ifl1 = 0

    ! Open L1 RAD files
    numFiles = 0
    do L1FileHandle = mlspcf_l1b_rad_start, mlspcf_l1b_rad_end
      returnStatus = InitializeMLSFile(L1BFile, content = 'l1brad', &
        & type=l_hdf, access=DFACC_RDONLY)
      L1BFile%PCFIDRange%Top = L1FileHandle
      L1BFile%PCFIDRange%Bottom = L1FileHandle
      L1BFile%hdfVersion = WILDCARDHDFVERSION
      call mls_openFile(L1BFile, returnStatus)
      if ( returnStatus == 0 ) then
        numFiles = addFileToDatabase(filedatabase, L1BFile)
        ifl1 = ifl1 + 1
      end if

    end do ! L1FileHandle = mlspcf_l1b_rad_start, mlspcf_l1b_rad_end

    if ( ifl1 == 0 .AND. TOOLKIT ) &
      &  call announce_error ( "Could not find any L1BRAD files" )

    ! Open L1OA File
    
    returnStatus = InitializeMLSFile(L1BFile, content = 'l1boa', &
      & type=l_hdf, access=DFACC_RDONLY)
    L1BFile%PCFIDRange%Top = mlspcf_l1b_oa_start
    L1BFile%PCFIDRange%Bottom = mlspcf_l1b_oa_start
    call mls_openFile(L1BFile, returnStatus)
    if ( returnStatus == PGS_S_SUCCESS ) then

        if( switchDetail(switches, 'pro') > -1 ) then  
          call announce_success(L1BFile%name, 'l1boa', &                     
          & hdfVersion=L1BFile%HDFVersion)                    
        end if
        numFiles = AddFileToDatabase(filedatabase, L1BFile)

        ! read l1boa attributes

        l1bFlag = 0
        call ReadL1BAttribute(l1bFile, OrbNum, 'OrbitNumber', &
          & l1bFlag)
        if (l1bFlag == -1) then
           GlobalAttributes%OrbNum = -1
        else
           GlobalAttributes%OrbNum = OrbNum
        end if
        call ReadL1BAttribute(l1bFile, OrbPeriod, 'OrbitPeriod', &
          & l1bFlag)
        if (l1bFlag == -1) then
           GlobalAttributes%OrbPeriod = -1.0
        else
           GlobalAttributes%OrbPeriod = OrbPeriod
        end if

    else if ( TOOLKIT ) then
      call announce_error ( "Could not find L1BOA file" )
    end if

    ! We'll possibly need the first and last MAF counter numbers, especially
    ! if gaps occur
    ! For now, just look for them in l1boa
    ! Later you may look also in l1brad files
    GlobalAttributes%LastMAFCtr = FindMaxMAF ( (/l1BFile/), &
      & GlobalAttributes%FirstMAFCtr )
    ! Get the Start and End Times from PCF

    returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_CCSDSStartId, &
                                           CCSDSStartTime)
    if ( returnstatus /= PGS_S_SUCCESS .and. TOOLKIT ) then
      call announce_error ( "Missing pcf param: CCSDSStartTime" )
    end if

    returnStatus = pgs_td_utctotai (CCSDSStartTime, processingrange%starttime)
    !   ??? Is PGSTD_E_NO_LEAP_SECS an OK status ???
    if ( returnstatus /= PGS_S_SUCCESS .and. &
      &  returnstatus /= PGSTD_E_NO_LEAP_SECS ) &
        & call announce_error ( "Could not convert UTC Start time to TAI" )
    if ( returnstatus == PGSTD_E_NO_LEAP_SECS ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'No leap second information' )

   returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_CCSDSEndId, &
                                          CCSDSEndTime)
    if ( returnstatus /= PGS_S_SUCCESS .and. TOOLKIT ) then
      call announce_error ( "Missing pcf param: CCSDSEndTime" )
    end if

    returnStatus = pgs_td_utctotai (CCSDSEndTime, processingrange%endtime)
    !   ??? Is PGSTD_E_NO_LEAP_SECS an OK status ???
    if ( returnstatus /= PGS_S_SUCCESS .and. &
      & returnstatus /= PGSTD_E_NO_LEAP_SECS) &
        & call announce_error ( "Could not convert UTC End time to TAI" )

    l2pcf%startutc = CCSDSStartTime
    l2pcf%endutc = CCSDSEndTime

    ! Here's where we define the non-time components of l2pcf

    returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_PGEVersion, &
                                          l2pcf%PGEVersion)
    if ( returnstatus /= PGS_S_SUCCESS ) then
      call announce_error ( "Missing pcf param: output version" )
    end if

    returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_Cycle, RunID ) ! l2pcf%cycle)
    if ( len_trim(RunID) < 5 ) then
      l2pcf%cycle = RunID
    else
      l2pcf%cycle = 'c01'
      l2pcf%RunID = RunID
    endif

    if ( returnstatus /= PGS_S_SUCCESS ) then
      call announce_error ( "Missing pcf param: cycle" )
    end if

    returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_col_spec_keys, col_species_keys)

    if ( returnstatus /= PGS_S_SUCCESS ) then
      call announce_error ( "Missing pcf param: col_spec_keys", &
        & forgiveable=.true. )
    end if

    returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_col_spec_hash, col_species_hash)

    if ( returnstatus /= PGS_S_SUCCESS ) then
      call announce_error ( "Missing pcf param: col_spec_hash", &
        & forgiveable=.true. )
    end if

    returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_spec_keys, l2pcf%spec_keys)

    if ( returnstatus /= PGS_S_SUCCESS ) then
      call announce_error ( "Missing pcf param: spec_keys" )
    end if

    returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_spec_hash, l2pcf%spec_hash)

    if ( returnstatus /= PGS_S_SUCCESS ) then
      call announce_error ( "Missing pcf param: spec_hash" )
    end if

    if ( .NOT. MCFCASESENSITIVE ) then
      l2pcf%spec_hash = LowerCase(l2pcf%spec_hash)
      l2pcf%spec_keys = LowerCase(l2pcf%spec_keys)
    end if

   ! If using PCF, it will be difficult to set switches from command-line
   ! so the following allows us to do so inside the PCF
   ! This way we can still dump vectors, monitor iterations or activities
   ! But don't do this if you already set them on the command line
    if ( ALWAYSGETSWITCHESCONFIG .or. switches == ' ' ) then
      returnStatus = pgs_pc_getConfigData(mlspcf_l2_param_switches, &
        & extra_switches)
      if ( returnstatus == PGS_S_SUCCESS .and. extra_switches /= ' ' ) &
        & switches = catLists(trim(switches), trim(extra_switches))
    end if

! Get the name of the log file from the PCF

    version = 1

    returnStatus = Pgs_pc_getReference(MLSPCF_LOG, version, name)
    if ( returnStatus /= PGS_S_SUCCESS ) then
      call announce_error ( "Error retrieving log file name from PCF" )
    elseif ( returnStatus == PGS_S_SUCCESS) then
      indx = INDEX(name, '/', .TRUE.)
      l2pcf%logGranID = name(indx+1:)
    else
      l2pcf%logGranID = '(not found)'
    end if
 
    ! Store appropriate user input as global attributes
    GlobalAttributes%StartUTC = l2pcf%StartUTC
    GlobalAttributes%EndUTC = l2pcf%EndUTC
    GlobalAttributes%PGEVersion = l2pcf%PGEVersion
    call utc_to_yyyymmdd(GlobalAttributes%StartUTC, returnStatus, &
      & GlobalAttributes%GranuleYear, GlobalAttributes%GranuleMonth, &
      & GlobalAttributes%GranuleDay) 
    call FillTAI93Attribute
    
    ! Get name of Parallel Staging file (not any more)
    call MLSMessage ( MLSMSG_Warning, ModuleName, &
      & 'This version does not use staging file for slave Join commands' )
    ! version = 1

    if ( error /= 0 ) &
      & call MLSMessage(MLSMSG_Error,ModuleName, &
        & 'Problem with open_init section')

   Details = switchDetail(switches, 'pcf') - 2 ! -3 means don't dump
   if ( levels(gen) > 0 .or. details > -3 .or. BeVerbose('pcf', -1) ) &
      & call Dump_open_init ( filedatabase, &
          & CCSDSEndTime, CCSDSStartTime, processingrange, details )
   if ( timing ) call sayTime         
   call trace_end ( "OpenAndInitialize", cond=toggle(gen) )

  contains
    subroutine SayTime
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      end if
      call output ( "Timing for open_init = " )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime
  end subroutine OpenAndInitialize

! =====     Private Procedures     =====================================

  ! ---------------------------------------------  announce_success  -----
  subroutine announce_success ( Name, l1_type, hdfVersion )
    character(LEN=*), intent(in)   :: Name
    character(LEN=*), intent(in)   :: l1_type
    integer, optional,  intent(in) :: hdfVersion

    call output ( 'Level 1 product type : ' )
    call output ( trim(l1_type), advance='no')
    if ( present(hdfVersion) ) then
      call blanks(4)
      call output ( 'hdf ' )
      call output ( hdfVersion, advance='yes')
    else
      call output ( ' ', advance='yes')
    end if
    call blanks(15)
    call output ( 'name : ' )
    call blanks(8)
    call output ( trim(Name), advance='yes')
  end subroutine announce_success

  ! ------------------------------------------  Dump_open_init  -----
  subroutine Dump_open_init ( filedatabase, &
    & CCSDSEndTime, CCSDSStartTime, processingrange, Details )
  
    ! Dump info obtained during OpenAndInitialize:
    ! L1B databse
    ! L1OA file
    ! Start and end times
    ! output version
    ! cycle number
    ! logfile name
  
    use L1BDATA, only: L1BDATA_T, NAME_LEN, &
      & ASSEMBLEL1BQTYNAME, DEALLOCATEL1BDATA, DUMP, READL1BDATA
    use L2AUXDATA, only: MAXSDNAMESBUFSIZE
    use L2GPDATA, only: COL_SPECIES_KEYS, COL_SPECIES_HASH
    use MLSFILES, only: HDFVERSION_4       
    use MLSHDF5, only: GETALLHDF5DSNAMES
    use OUTPUT_M, only: REVERTOUTPUT, SWITCHOUTPUT
    use WRITEMETADATA, only: L2PCF

    ! Arguments
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    character(len=CCSDSlen) CCSDSEndTime
    character(len=CCSDSlen) CCSDSStartTime
    type (TAI93_Range_T) :: processingRange ! Data processing range
    integer, intent(in) :: Details

    ! Local

    integer :: i, NoMAFs, IERR
    character (len=*), parameter :: time_format='(1pD18.12)'
    character (len=*), parameter :: l1b_quant_name='R1A:118.B1F:PT.S0.FB25-1'
    type (L1BData_T) :: l1bDataSet   ! L1B dataset

    integer ::  hdfVersion
    integer :: ifl1
    integer :: item
    character(len=Name_Len) :: l1bItemName
    character(len=MAXSDNAMESBUFSIZE) :: l1bitemlist
    integer :: nL1bItems
    logical, parameter :: countEmpty = .true.
    type (MLSFile_T), pointer ::     L1BFile

    ! Begin                                                                       
    call output ( '============ Open Initialize ============', advance='yes' )
    call output ( '        (Data entered via pcf)', advance='yes' )
    call output ( ' ', advance='yes' )
    call output ( 'L1B database:', advance='yes' )
  
   if(associated(filedatabase)) then
     if ( specialDumpFile /= ' ' .and. Details > -2 ) &
       & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
    ifl1 = 0
    do i = 1, size(filedatabase)
      L1BFile => filedatabase(i)
      if ( L1BFile%content == 'l1brad' ) then
         call output ( 'PCF id:   ' )                                       
         call output ( L1BFile%PCFID, advance='yes' )
         call output ( 'fileid:   ' )
         call output ( L1BFile%FileID%f_id, advance='yes' )
         call output ( 'name:   ' )                                       
         call output ( TRIM(L1BFile%name), advance='yes' )
         call output ( 'hdf version:   ' )                                       
         call output ( L1BFile%hdfVersion, advance='yes' )
         if ( Details > -2 .and. hdfVersion == HDFVERSION_4 ) then
           l1bItemName = AssembleL1BQtyName ( l1b_quant_name, hdfVersion, .false. )
           call ReadL1BData ( L1BFile, l1bItemname, L1bDataSet, &
            & NoMAFs, IERR, NeverFail=.true. )
           if ( IERR == 0 ) then
             call Dump ( l1bDataSet, Details )
             call DeallocateL1BData ( l1bDataSet )
           else
             call output ( 'Error number  ' )
             call output ( IERR )
             call output ( ' while reading quantity named  ' )
             call output ( trim(l1b_quant_name) )
           end if
         elseif( Details > -2 ) then
           call GetAllHDF5DSNames (L1BFile%FileID%f_id, '/', l1bitemlist)
           nL1bItems = NumStringElements(trim(l1bitemlist), countEmpty)
           do item=1, nL1bItems
             call GetStringElement (trim(l1bitemlist), l1bItemName, item, &
               & countEmpty )
             if ( index('/PCF,/LCF', trim(l1bItemName)) > 0 ) cycle
             call ReadL1BData ( L1BFile, l1bItemName, L1bDataSet, &
              & NoMAFs, IERR, NeverFail=.true. )
             if ( IERR == 0 ) then
               call Dump ( l1bDataSet, Details )
               call DeallocateL1BData ( l1bDataSet )
             else
               call output ( 'Error number  ' )
               call output ( IERR )
               call output ( ' while reading quantity named  ' )
               call output ( trim(l1bItemName) )
             end if
             ! if ( item > 2 ) stop
           enddo
        end if
      elseif ( L1BFile%content == 'l1boa' ) then
        call output ( 'L1OA file:', advance='yes' )
  
        call output ( 'PCFid:   ' )
        call output ( L1BFile%PCFId, advance='yes' )
        call output ( 'fileid:   ' )
        call output ( L1BFile%FileID%f_id, advance='yes' )
        call output ( 'name:   ' )
        call output ( TRIM(L1BFile%Name), advance='yes' )
      end if
    enddo

     if ( specialDumpFile /= ' ' .and. Details > -2 ) &
       & call revertOutput

   else
    call output ( '(null database)', advance='yes' )

   end if

    call output ( 'Start Time:   ' )
    call output ( CCSDSStartTime, advance='yes' )

    call output ( 'End Time:     ' )
    call output ( CCSDSEndTime, advance='yes' )

    call output ( 'Start Time (tai):   ' )
    call output ( processingrange%starttime, format=time_format, advance='yes' )

    call output ( 'End Time (tai):     ' )
    call output ( processingrange%endtime, format=time_format, advance='yes' )

    call output ( 'Processing Range:     ' )
    call output ( processingrange%endtime-processingrange%starttime, advance='yes' )

    call output ( 'PGE version:   ' )
    call output ( l2pcf%PGEVersion, advance='yes' )

    call output ( 'cycle:   ' )
    call output ( l2pcf%cycle, advance='yes' )

    call output ( 'RunID:   ' )
    call output ( l2pcf%RunID, advance='yes' )

    call output ( 'Log file name:   ' )
    call output ( TRIM(l2pcf%logGranID), advance='yes' )

    call output ( 'l2gp column species name keys:   ' )
    call output ( TRIM(col_species_keys), advance='yes' )

    call output ( 'corresponding units hash:   ' )
    call output ( TRIM(col_species_hash), advance='yes' )

    call output ( 'l2gp species name keys:   ' )
    call output ( TRIM(l2pcf%spec_keys), advance='yes' )

    call output ( 'corresponding mcf hash:   ' )
    call output ( TRIM(l2pcf%spec_hash), advance='yes' )

  end subroutine Dump_open_init

  ! ---------------------------------------------  Announce_Error  -----
  subroutine Announce_Error ( Full_message, Use_toolkit, &
    & Error_number, forgiveable )
    ! Arguments

    ! integer, intent(in) :: Lcf_where
    character(LEN=*), intent(in) :: Full_message
    logical, intent(in), optional :: Use_toolkit
    integer, intent(in), optional :: Error_number
    logical, intent(in), optional :: forgiveable

    ! Local
    logical :: Just_print_it
    logical, parameter :: Default_output_by_toolkit = .true.

    just_print_it = .not. default_output_by_toolkit
    if ( present(use_toolkit) ) just_print_it = .not. use_toolkit
    if ( present(forgiveable) ) just_print_it = just_print_it .or. forgiveable

    if ( .not. just_print_it ) then
      error = max(error,1)
      call output ( " The following error occurred:", advance='yes', &
       & from_where=ModuleName)
      call output ( trim(full_message), advance='yes', &
        & from_where=ModuleName)
      if ( present(error_number) ) then
        call output ( 'Error number ', advance='no' )
        call output ( error_number, places=9, advance='yes' )
      end if
    else
      call output ( '*** ' )
      if ( present(forgiveable) ) then
        call output ( 'Warning in module' )
      else
        call output ( 'Error in module ' )
      endif
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

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Open_Init

!=============================================================================

!
! $Log$
! Revision 2.105  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.104  2013/12/05 01:43:48  pwagner
! Read RunID component from cycle field of pcf; started using BeVerbose
!
! Revision 2.103  2013/08/30 02:45:50  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.102  2013/08/17 02:55:23  vsnyder
! Regularized trace usage
!
! Revision 2.101  2012/08/16 17:49:40  pwagner
! Simplify; remove unused stuff
!
! Revision 2.100  2012/04/05 20:14:39  pwagner
! Erased all traces of parallel staging file
!
! Revision 2.99  2012/03/28 20:09:38  pwagner
! Issue warning--slave tasks lost ability to join quantities
!
! Revision 2.98  2011/06/29 21:49:43  pwagner
! Some cases may safely omit l1b files
!
! Revision 2.97  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.96  2007/10/04 20:43:12  vsnyder
! Remove unused symbols
!
! Revision 2.95  2007/06/21 00:54:08  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.94  2006/09/21 18:50:33  pwagner
! Uncertain why or whether this was necessary but did it anyway
!
! Revision 2.93  2006/07/07 23:11:19  pwagner
! Removed already-commented-out souvenirs
!
! Revision 2.92  2006/03/15 23:52:24  pwagner
! Removed InputVersion component from PCF, l2cf
!
! Revision 2.91  2006/02/16 00:16:01  pwagner
! switchDetail instead of index
!
! Revision 2.90  2006/02/10 21:16:20  pwagner
! dumps may go to special dumpfile
!
! Revision 2.89  2006/01/27 01:03:29  pwagner
! Changed wording when merely a warning
!
! Revision 2.88  2006/01/19 00:31:54  pwagner
! reads col_spec_keys and _hash from PCF
!
! Revision 2.87  2005/09/22 23:38:12  pwagner
! date conversion procedures and functions all moved into dates module
!
! Revision 2.86  2005/07/21 23:45:21  pwagner
! Removed unused l1b fileinfo fields from l2pcf
!
! Revision 2.85  2005/07/12 17:17:46  pwagner
! Dropped global attribute InputVersion
!
! Revision 2.84  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.83  2005/06/14 20:39:25  pwagner
! Many changes to accommodate the new fields in MLSFile_T
!
! Revision 2.82  2005/05/31 17:51:17  pwagner
! Began switch from passing file handles to passing MLSFiles
!
! Revision 2.81  2004/12/28 00:23:03  vsnyder
! Remove unreferenced use names
!
! Revision 2.80  2004/12/15 23:35:59  pwagner
! Should not warn about missing switches setting in PCF
!
! Revision 2.79  2004/08/16 17:14:42  pwagner
! Obtains First,LastMAFCtr global attributes from l1boa file
!
! Revision 2.78  2004/08/04 23:19:58  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.77  2004/06/29 00:10:17  pwagner
! Exploit catlist function
!
! Revision 2.76  2003/10/29 00:07:17  pwagner
! GlobalAttributes%ProcessLevel assigned appropriate value
!
! Revision 2.75  2003/10/09 23:41:06  pwagner
! Ignore illegal PCFids; extra_switches ,-separated
!
! Revision 2.74  2003/09/12 16:28:56  cvuu
! Read OrbitNumber and OrbitPeriod attributes from L1BOA
!
! Revision 2.73  2003/07/07 23:50:05  pwagner
! Now uses saved variable L2pcf from writeMetaData
!
! Revision 2.72  2003/06/09 22:49:34  pwagner
! Reduced everything (PCF, PUNISH.., etc.) to TOOLKIT
!
! Revision 2.71  2003/05/29 17:55:23  pwagner
! Can get name of parallel%stagingFile from PCF
!
! Revision 2.70  2003/02/27 21:55:14  pwagner
! Calls FillTAI93Attribute
!
! Revision 2.69  2003/02/01 00:38:33  pwagner
! Defines Global Attributes from user data
!
! Revision 2.68  2002/12/11 22:18:48  pwagner
! broadened error check on sd_id to any value lt 1
!
! Revision 2.67  2002/11/13 01:10:40  pwagner
! Actually reads hdf5 radiances
!
! Revision 2.66  2002/10/08 17:36:22  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.65  2002/10/03 23:02:48  pwagner
! Now opens l1b files with mls_io_gen_openF
!
! Revision 2.64  2002/09/27 23:47:23  pwagner
! Now calls mls_io_gen_closeF rather than mls_sfend
!
! Revision 2.63  2002/08/28 22:28:37  pwagner
! Saves PCFids for l1boa, l1brad for input pointer metadata
!
! Revision 2.62  2002/08/21 02:23:39  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.61  2002/01/26 00:09:14  pwagner
! Convert hdf routines to mlsfiles equivalents
!
! Revision 2.60  2002/01/24 00:14:26  pwagner
! Proclaims level 1 files as input; comments concerning hdf5 conversion
!
! Revision 2.59  2001/12/16 00:57:46  livesey
! Add warning about leap seconds
!
! Revision 2.58  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.57  2001/11/01 21:06:17  pwagner
! Fixed if block; added toc
!
! Revision 2.56  2001/10/30 00:33:24  pwagner
! Renamed Dump_L1B_database Dump_open_init; switches now checked for pcf[n]
!
! Revision 2.55  2001/10/25 23:33:59  pwagner
! Initializes InputVersion, PGEVersion to blank
!
! Revision 2.54  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.53  2001/08/09 00:23:05  pwagner
! pcf param for extra switches added; read
!
! Revision 2.52  2001/05/24 20:34:50  pwagner
! More forgiving of faulty pcf
!
! Revision 2.51  2001/05/17 22:31:54  pwagner
! Simpler info instead of bogus dump if no pcf
!
! Revision 2.50  2001/05/17 00:27:55  pwagner
! Better defaults, dumps
!
! Revision 2.49  2001/05/11 23:43:35  pwagner
! Added some initializations
!
! Revision 2.48  2001/05/10 18:22:57  pwagner
! Improved Dump_L1B_database
!
! Revision 2.47  2001/05/09 23:32:15  pwagner
! Gets toolkit funs from SDPToolkit
!
! Revision 2.46  2001/05/07 23:28:19  pwagner
! Fixed dump..
!
! Revision 2.45  2001/05/06 20:55:25  pwagner
! [De]Allocates l1binfo%filenames along with ids
!
! Revision 2.44  2001/05/04 22:56:40  pwagner
! Detachable from Toolkit
!
! Revision 2.43  2001/05/04 17:09:03  pwagner
! Get MAXNUML1BRADIDS, ILLEGALL1BRADID from global_settings
!
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

