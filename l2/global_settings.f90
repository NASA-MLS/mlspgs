! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GLOBAL_SETTINGS

  use EXPR_M, only: EXPR   
  use ForwardModelConfig, only: AddForwardModelConfigToDatabase, &
    & ForwardModelConfig_T
  use ForwardModelSupport, only: ConstructForwardModelConfig, &
    & ForwardModelGlobalSetup
  use INIT_TABLES_MODULE, only: L_TRUE, P_ALLOW_CLIMATOLOGY_OVERLOADS, &
    & P_INPUT_VERSION_STRING, P_OUTPUT_VERSION_STRING, P_VERSION_COMMENT, &
    & S_FORWARDMODEL, S_ForwardModelGlobal, S_TIME, S_VGRID, F_FILE, &
    & P_CYCLE, P_STARTTIME, P_ENDTIME, &
    & S_L1BRAD, S_L1BOA
  use L1BData, only: l1bradSetup, l1boaSetup, ReadL1BData, L1BData_T, NAME_LEN, &
    & DeallocateL1BData
  use L2GPData, only: L2GPDATA_T
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: R8, NameLen, L1BInfo_T, TAI93_Range_T, FileNameLen
  use MLSL2Options, only: PCF
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate
  use MLSPCF2, only: MLSPCF_L1B_RAD_END, MLSPCF_L1B_RAD_START
  use MLSStrings, only: unquote, hhmmss_value
  use MoreTree, only: GET_FIELD_ID, GET_SPEC_ID
  use Output_m, only: Output
  use String_Table, only: Get_String
  use TOGGLES, only: GEN, LEVELS, SWITCHES, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, SUB_ROSA, SUBTREE, &
  & DUMP_TREE_NODE, SOURCE_REF
  use TREE_TYPES, only: N_EQUAL, N_NAMED
  use VGrid, only: CreateVGridFromMLSCFInfo
  use VGridsDatabase, only: AddVGridToDatabase, Dump, VGrid_T
  use WriteMetadata, only: PCFData_T

  implicit NONE

  private

  public :: SET_GLOBAL_SETTINGS

  integer, public, parameter :: ILLEGALL1BRADID=-1      ! something sfstart should catch
  integer, public, parameter :: MAXNUML1BRADIDS=&
  & mlspcf_l1b_rad_end-mlspcf_l1b_rad_start+1   ! In case more than one

  logical, public :: ALLOW_CLIMATOLOGY_OVERLOADS = .false.
  integer, public :: INPUT_VERSION_STRING = 0     ! Sub_rosa index
  integer, public :: OUTPUT_VERSION_STRING = 0    ! Sub_rosa index
  integer, public :: VERSION_COMMENT = 0          ! Sub_rosa index

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  integer, private :: ERROR

contains

  subroutine SET_GLOBAL_SETTINGS ( ROOT, ForwardModelConfigDatabase, &
    & VGrids, l2gpDatabase, l2pcf, processingRange, l1bInfo )

    integer, intent(in) :: ROOT    ! Index of N_CF node in abstract syntax tree
    type(ForwardModelConfig_T), dimension(:), pointer :: &
      & ForwardModelConfigDatabase
    type ( vGrid_T ), pointer, dimension(:) :: VGrids
    type ( l2gpData_T), dimension(:), pointer :: L2GPDATABASE
    type (TAI93_Range_T) :: processingRange ! Data processing range
    type (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset
    type (L1BData_T) :: l1bField ! L1B data
    type(PCFData_T) :: l2pcf
    
    logical :: GOT(2) = .false.
    integer :: I         ! Index of son of root
    integer :: KEY       ! A P_... parameter from Init_Tables_Module
    integer :: L1BFLAG
    real(r8) :: MINTIME, MAXTIME        ! Time Span in L1B file data
    integer :: NAME      ! Sub-rosa index of name of vGrid or hGrid
    integer :: NOMAFS             ! Number of MAFs of L1B data read
    integer :: returnStatus             ! non-zero means trouble
    integer :: SON       ! Son of root
    integer :: sub_rosa_index
    logical :: TIMING    ! For S_Time
    real :: T1, T2       ! For S_Time
    integer :: UNITS(2)  ! Units of expression
    real(r8) :: VALUE(2)   ! Value of expression
    real(r8) :: start_time_from_1stMAF, end_time_from_1stMAF

    character(LEN=NameLen) :: name_string
    character (len=name_len) :: QUANTITY
    character(LEN=*), parameter :: time_conversion='(F32.0)'

    timing = .false.
    
   error = 0

    if ( toggle(gen) ) call trace_begin ( 'SET_GLOBAL_SETTINGS', root )

    do i = 2, nsons(root)-1 ! Skip names at beginning and end of section
      son = subtree(i,root)
      if ( node_id(son) == n_equal ) then
        sub_rosa_index = sub_rosa(subtree(2,son))
        select case ( decoration(subtree(1,son)) )
        case ( p_allow_climatology_overloads )
          allow_climatology_overloads = decoration(subtree(2,son)) == l_true
        case ( p_input_version_string )
          input_version_string = sub_rosa_index
          call get_string ( input_version_string, l2pcf%inputVersion, strip=.true. )
          if ( pcf ) &
            & call announce_error(0, &
            & '*** l2cf overrides pcf for input version ***', &
            & just_a_warning = .true.)
        case ( p_output_version_string )
          output_version_string = sub_rosa_index
          call get_string ( output_version_string, l2pcf%PGEVersion, strip=.true. )
          if ( pcf ) &
            & call announce_error(0, &
            & '*** l2cf overrides pcf for PGE version ***', &
            & just_a_warning = .true.)
        case ( p_version_comment )
          version_comment = sub_rosa_index
        case ( p_cycle )
          call get_string ( sub_rosa_index, l2pcf%cycle, strip=.true. )
          if ( pcf ) &
            & call announce_error(0, &
            & '*** l2cf overrides pcf for cycle number ***', &
            & just_a_warning = .true.)
!        case ( p_ccsdsstarttime )
!          call get_string ( sub_rosa_index, l2pcf%startutc, strip=.true. )
!        case ( p_ccsdsendtime )
!          call get_string ( sub_rosa_index, l2pcf%endutc, strip=.true. )
        case ( p_starttime )
          got(1) = .true.
          call get_string ( sub_rosa_index, name_string, strip=.true. )
          if ( index(name_string, ':') > 0 ) then
            start_time_from_1stMAF = hhmmss_value(name_string, error)
          else
            call expr ( subtree(2,son), units, value )
            start_time_from_1stMAF = value(1)
          endif
          if ( pcf ) &
            & call announce_error(0, &
            & '*** l2cf overrides pcf for start time ***', &
            & just_a_warning = .true.)
 !         print *, 'starttime'
 !         print *, trim(name_string)
 !         print *, trim(unquote(name_string))
 !         name_string = unquote(name_string)
 !         read(name_string, time_conversion) &
 !         &           processingrange%starttime
        case ( p_endtime )
          got(2) = .true.
          call get_string ( sub_rosa_index, name_string, strip=.true. )
          if ( index(name_string, ':') > 0 ) then
            end_time_from_1stMAF = hhmmss_value(name_string, error)
          else
            call expr ( subtree(2,son), units, value )
            end_time_from_1stMAF = value(1)
          endif
          if ( pcf ) &
            & call announce_error(0, &
            & '*** l2cf overrides pcf for end time ***', &
            & just_a_warning = .true.)
 !         call get_string ( sub_rosa_index, name_string, strip=.true. )
 !         print *, 'endtime'
 !         print *, trim(name_string)
 !         print *, trim(unquote(name_string))
 !         name_string = unquote(name_string)
 !         read(name_string, time_conversion) &
 !         &           processingrange%endtime
        case default
         call announce_error(son, 'unrecognized global settings parameter')
        end select
      else
        if ( node_id(son) == n_named ) then
          name = sub_rosa(subtree(1,son))
          son = subtree(2,son)
        else
          name = 0
        end if
        select case ( get_spec_id(son) )
        case ( s_forwardModelGlobal ) !??? Begin temporary stuff for l2load
          call forwardModelGlobalSetup ( son, returnStatus )
          error = max(error, returnStatus)
        case ( s_forwardModel )
          call decorate (son, AddForwardModelConfigToDatabase ( &
            & forwardModelConfigDatabase, ConstructForwardModelConfig ( son, vGrids ) ) )
        case ( s_vgrid )
          call decorate ( son, AddVGridToDatabase ( vGrids, &
            & CreateVGridFromMLSCFInfo ( name, son, l2gpDatabase ) ) )
        case ( s_l1brad )
          call l1bradSetup ( son, l1bInfo, F_FILE, &
          & MAXNUML1BRADIDS, ILLEGALL1BRADID )
          if ( pcf ) &
            & call announce_error(0, &
            & '*** l2cf overrides pcf for L1B Rad. file(s) ***', &
            & just_a_warning = .true.)
        case ( s_l1boa )
          call l1boaSetup ( son, l1bInfo, F_FILE )
          if ( pcf ) &
            & call announce_error(0, &
            & '*** l2cf overrides pcf for L1BOA file ***', &
            & just_a_warning = .true.)

        case ( s_time )
          if ( timing ) then
            call sayTime
          else
            call cpu_time ( t1 )
            timing = .true.
          end if
        case default
         call announce_error(son, 'unrecognized global settings spec')
        end select
      end if
    end do

   ! add maf offsets to start, end times
   ! This is optional way to define processingRange if using PCF
   ! It becomes mandatory if not using PCF
   if(got(1) .or. got(2) .or. .not. PCF) then
   
   ! 1st--check that have L1BOA
     if(l1bInfo%L1BOAID == ILLEGALL1BRADID) then
       call announce_error(son, &
       & 'L1BOA file required by global data--but not set')
     endif
     quantity = 'MAFStartTimeTAI'
     call ReadL1BData ( l1bInfo%l1boaID, quantity, l1bField, noMAFs, &
          & l1bFlag)
      if ( l1bFlag==-1) then
            call announce_error(son, &
          & 'unrecognized MAFStarttimeTAI in L1BOA file')
           minTime = 0.
           maxTime = 0.
      else
           minTime = l1bField%dpField(1,1,1)
           maxTime = l1bField%dpField(1,1,noMAFs) ! This is start time of last MAF
      endif
      call DeallocateL1BData ( l1bField )
   endif
   if(got(1)) then
      processingrange%starttime = minTime + start_time_from_1stMAF
   elseif(.not. PCF) then
      processingrange%starttime = minTime
   endif

   if(got(2)) then
      processingrange%endtime = minTime + end_time_from_1stMAF
   elseif(.not. PCF) then
      processingrange%endtime = maxTime + 1.0
   endif

   if( levels(gen) > 0 .or. &
   & index(switches, 'glo') /= 0 ) &
   & call dump_global_settings( l2pcf, processingRange, l1bInfo )

    if ( error /= 0 ) &
      & call MLSMessage(MLSMSG_Error,ModuleName, &
        & 'Problem with global settings section')

    if ( toggle(gen) ) then
      if (  levels(gen) > 0 .or. index(switches, 'V') /= 0 ) &
        & call dump ( vgrids, details=levels(gen)-1+min(index(switches, 'V'),1) )
      call trace_end ( 'SET_GLOBAL_SETTINGS' )
    end if
    if ( timing ) call sayTime

  contains

    ! --------------------------------------------------  SayTime  -----
    subroutine SayTime
      call cpu_time ( t2 )
      call output ( "Timing for GlobalSettings = " )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

  end subroutine SET_GLOBAL_SETTINGS

  ! ------------------------------------------  dump_global_settings  -----
  subroutine dump_global_settings ( l2pcf, processingRange, l1bInfo )
  
    ! Dump info obtained during OpenAndInitialize and global_settings:
    ! L1B databse
    ! L1OA file
    ! Start and end times
    ! output version
    ! cycle number
    ! logfile name
  
    integer :: num_l1b_files = MAXNUML1BRADIDS

    ! Arguments
    type (L1BInfo_T) :: l1bInfo   ! File handles etc. for L1B dataset
    type(PCFData_T) :: l2pcf
    type (TAI93_Range_T) :: processingRange ! Data processing range
!    character(len=CCSDSlen) CCSDSEndTime
!    character(len=CCSDSlen) CCSDSStartTime
	
    ! Local

    character (LEN=FileNameLen) :: physicalFilename
    integer :: i, returnStatus, version
    character (len=*), parameter :: time_format='(1pD18.12)'

    ! Begin
    version = 1

    call output ( '============ Global Settings ============', advance='yes' )
    call output ( ' ', advance='yes' )
  
    call output ( 'L1B database:', advance='yes' )
  
   if(associated(l1bInfo%L1BRADIDs)) then
    if ( num_l1b_files > 0 ) then
      do i = 1, num_l1b_files
      if(l1bInfo%L1BRADIDs(i) /= ILLEGALL1BRADID) then
  	    call output ( 'fileid:   ' )
	    call output ( l1bInfo%L1BRADIDs(i), advance='yes' )
      	call output ( 'name:   ' )
    	   call output ( TRIM(l1bInfo%L1BRADFileNames(i)), advance='yes' )
      endif
      end do

    else
      call output ( '(empty database)', advance='yes' )
    end if

   else
    call output ( '(null database)', advance='yes' )

   endif

    call output ( ' ', advance='yes' )
    call output ( 'L1OA file:', advance='yes' )
  
      if(l1bInfo%L1BOAID /= ILLEGALL1BRADID) then
      call output ( 'fileid:   ' )
      call output ( l1bInfo%L1BOAID, advance='yes' )
      call output ( 'name:   ' )
      call output ( TRIM(l1bInfo%L1BOAFileName), advance='yes' )
    else
      call output ( '(file unknown)', advance='yes' )
    end if

    call output ( ' ', advance='yes' )
    call output ( 'Start Time:   ' )
    call output ( l2pcf%startutc, advance='yes' )

    call output ( 'End Time:     ' )
    call output ( l2pcf%endutc, advance='yes' )

    call output ( 'Start Time (tai):   ' )
    call output ( processingrange%starttime, format=time_format, advance='yes' )

    call output ( 'End Time (tai):     ' )
    call output ( processingrange%endtime, format=time_format, advance='yes' )

    call output ( 'Processing Range:     ' )
    call output ( processingrange%endtime-processingrange%starttime, advance='yes' )

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

    call output ( 'Allow climatology overloads?:   ' )
    call output ( allow_climatology_overloads, advance='yes' )

    call output ( ' ', advance='yes' )
    call output ( '============ End Global Settings ============', advance='yes' )

  end subroutine dump_global_settings

  ! ---------------------------------------------  Announce_Error  -----
  subroutine Announce_Error ( Lcf_where, Full_message, Use_toolkit, &
    & Error_number, just_a_warning )
  
    ! Arguments

    integer, intent(in) :: Lcf_where
    character(LEN=*), intent(in) :: Full_message
    logical, intent(in), optional :: Use_toolkit
    integer, intent(in), optional :: Error_number
    logical, intent(in), optional :: just_a_warning

    ! Local
    logical :: Just_print_it, my_warning
    logical, parameter :: Default_output_by_toolkit = .true.

    just_print_it = .not. default_output_by_toolkit
    if ( present(use_toolkit) ) just_print_it = .not. use_toolkit
    if ( present(just_a_warning) ) then
      my_warning = just_a_warning
    else
      my_warning = .false.
    endif

    if ( .not. just_print_it ) then
     if ( .not. my_warning ) then
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

      endif
      
      call output ( trim(full_message), advance='yes', &
        & from_where=ModuleName)

      if ( present(error_number) ) then
       if (my_warning) then
        call output ( 'Warning number ', advance='no' )
      else
        call output ( 'Error number ', advance='no' )
      endif

        call output ( error_number, places=9, advance='yes' )
      end if
    else

      if ( .not. my_warning) then
      call output ( '***Error in module ' )
      call output ( ModuleName, advance='yes' )
      endif
      
      call output ( trim(full_message), advance='yes' )
      if ( present(error_number) ) then
       if(my_warning) then
        call output ( 'Warning number ' )
      else
        call output ( 'Error number ' )
       endif
        call output ( error_number, advance='yes' )
      end if
    end if

!===========================
  end subroutine Announce_Error
!===========================

end module GLOBAL_SETTINGS

! $Log$
! Revision 2.39  2001/07/09 22:53:20  pwagner
! Obeys CloudForwardModel; for now same as ForwardModel
!
! Revision 2.38  2001/06/07 21:58:28  pwagner
! Added Copyright statement
!
! Revision 2.37  2001/05/30 23:56:23  livesey
! Changed for new L1BData
!
! Revision 2.36  2001/05/30 23:04:40  pwagner
! Gets returnStatus from forwardModelGlobalSetup
!
! Revision 2.35  2001/05/29 23:21:07  livesey
! Now uses ForwardModelSupport, not ForwardModelInterface
!
! Revision 2.34  2001/05/26 00:06:49  livesey
! Added call to DealloteL1BData
!
! Revision 2.33  2001/05/24 20:54:15  pwagner
! Deleted p_ccs..times
!
! Revision 2.32  2001/05/24 20:36:13  pwagner
! Warns if glob. stg. overrides pcf
!
! Revision 2.31  2001/05/17 00:29:03  pwagner
! Works without toolkit, PCF at last
!
! Revision 2.30  2001/05/15 23:46:32  pwagner
! Now optionally uses hhmmss_value
!
! Revision 2.29  2001/05/14 23:45:08  pwagner
! Start, end times now added to MAF offsets from L1BOA
!
! Revision 2.28  2001/05/11 23:44:43  pwagner
! Better dump; uses strip=TRUE
!
! Revision 2.27  2001/05/11 01:56:17  vsnyder
! Move the getting of sub_rosa_index
!
! Revision 2.26  2001/05/11 00:09:05  pwagner
! Gets p_.. from init_tables; unquotes strings
!
! Revision 2.25  2001/05/10 18:26:22  pwagner
! Improved dump_global_settings
!
! Revision 2.24  2001/05/09 23:35:04  pwagner
! Added dump_global_settings
!
! Revision 2.23  2001/05/04 17:15:36  pwagner
! Many added settings, esp. L1B files, so level2 can run w/o PCF
!
! Revision 2.22  2001/04/26 20:02:09  livesey
! Made l2pc database a saved array in L2PC_m
!
! Revision 2.21  2001/04/26 02:52:17  vsnyder
! Fix up CVS stuff
!
! Revision 2.20  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.19  2001/04/26 00:07:33  livesey
! Stuff to support reading of l2pc files
!
! Revision 2.18  2001/04/24 00:31:42  vsnyder
! Finish adding 'time' command
!
! Revision 2.17  2001/04/23 23:48:41  vsnyder
! Finish adding 'time' command
!
! Revision 2.16  2001/04/23 23:42:00  vsnyder
! Add 'time' command
!
! Revision 2.15  2001/04/21 01:25:54  livesey
! Now passes l2gpdatabase to more people who need it.
!
! Revision 2.14  2001/04/10 23:44:44  vsnyder
! Improve 'dump'
!
! Revision 2.13  2001/04/10 02:46:17  livesey
! Working version, no more FMI/TFMI
!
! Revision 2.12  2001/04/07 01:50:49  vsnyder
! Move some of VGrid to lib/VGridsDatabase.  Move ForwardModelConfig_T and
! some related stuff to fwdmdl/ForwardModelConfig.
!
! Revision 2.11  2001/03/28 22:00:33  livesey
! Interim version, now handles vGrids as part of forwardModelConfig
!
! Revision 2.10  2001/03/28 01:24:55  vsnyder
! Move vGrid from construct section to global settings section
!
! Revision 2.9  2001/03/17 03:24:23  vsnyder
! Work on forwardModelGlobalSetup
!
! Revision 2.8  2001/03/17 00:57:36  livesey
! Removed dump.
!
! Revision 2.7  2001/03/17 00:45:38  livesey
! Added ForwardModelConfigDatabase
!
! Revision 2.6  2001/03/09 02:30:13  vsnyder
! Allocate correct size for FMI and TFMI
!
! Revision 2.5  2001/03/09 00:24:30  vsnyder
! Do subscripts right
!
! Revision 2.4  2001/03/08 03:23:09  vsnyder
! More stuff to work with L2_Load
!
! Revision 2.3  2001/03/07 22:46:05  vsnyder
! Add temporary stuff for Zvi's "l2_load", which will wither away.
!
! Revision 2.2  2000/11/16 01:53:40  vsnyder
! Take timing back out.  Don't do it in sections that are only parameter settings.
!
! Revision 2.1  2000/11/16 01:45:25  vsnyder
! Implement timing.
!
! Revision 2.0  2000/09/05 18:57:05  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!
