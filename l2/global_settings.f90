module GLOBAL_SETTINGS

  use ForwardModelConfig, only: AddForwardModelConfigToDatabase, &
    & ForwardModelConfig_T
  use ForwardModelInterface, only: ConstructForwardModelConfig, &
    & ForwardModelGlobalSetup
  use HGrid, only: AddHGridToDatabase, CreateHGridFromMLSCFInfo, &
    & DestroyHGridDatabase, HGrid_T
  use INIT_TABLES_MODULE, only: L_TRUE, P_ALLOW_CLIMATOLOGY_OVERLOADS, &
    & P_INPUT_VERSION_STRING, P_OUTPUT_VERSION_STRING, P_VERSION_COMMENT, &
    & S_FORWARDMODEL, S_ForwardModelGlobal, S_TIME, S_VGRID, F_FILE!, &
!    & P_CYCLE, P_CCSDSSTARTTIME, P_CCSDSENDTIME, P_STARTTIME, P_ENDTIME, &
!    & S_L1BRAD, S_L1BOA
  use L1BData, only: l1bradSetup, l1boaSetup
  use L2GPData, only: L2GPDATA_T
  use MLSCommon, only: R8, NameLen, L1BInfo_T, TAI93_Range_T, FileNameLen
  use MLSL2Options, only: ECHO_GLOBAL_STNGS
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate
  use MLSPCF2, only: MLSPCF_L1B_RAD_END, MLSPCF_L1B_RAD_START
  use MoreTree, only: GET_FIELD_ID, GET_SPEC_ID
  use Output_m, only: Output
  use String_Table, only: Get_String
  use TOGGLES, only: GEN, LEVELS, SWITCHES, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, SUB_ROSA, SUBTREE
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
    type(PCFData_T) :: l2pcf
    
    integer :: I         ! Index of son of root
    integer :: NAME      ! Sub-rosa index of name of vGrid or hGrid
    integer :: SON       ! Son of root
    logical :: TIMING    ! For S_Time
    real :: T1, T2       ! For S_Time
    character(LEN=NameLen) :: name_string
! Just until init_tables_module is updated
   integer, parameter :: P_CYCLE=-99, P_CCSDSSTARTTIME=-98, &
   & P_CCSDSENDTIME=-97, &
   & P_STARTTIME=-96, P_ENDTIME=-95, &
    & S_L1BRAD=-94, S_L1BOA=-93

    timing = .false.
    
    if ( toggle(gen) ) call trace_begin ( 'SET_GLOBAL_SETTINGS', root )

    do i = 2, nsons(root)-1 ! Skip names at beginning and end of section
      son = subtree(i,root)
      if ( node_id(son) == n_equal ) then
        select case ( decoration(subtree(1,son)) )
        case ( p_allow_climatology_overloads )
          allow_climatology_overloads = decoration(subtree(2,son)) == l_true
        case ( p_input_version_string )
          input_version_string = sub_rosa(subtree(2,son))
        case ( p_output_version_string )
          output_version_string = sub_rosa(subtree(2,son))
        case ( p_version_comment )
          version_comment = sub_rosa(subtree(2,son))
        case ( p_cycle )
          call get_string ( sub_rosa(subtree(2,son)), l2pcf%cycle )
        case ( p_ccsdsstarttime )
          call get_string ( sub_rosa(subtree(2,son)), l2pcf%startutc )
        case ( p_ccsdsendtime )
          call get_string ( sub_rosa(subtree(2,son)), l2pcf%endutc )
        case ( p_starttime )
          call get_string ( sub_rosa(subtree(2,son)), name_string )
          read(name_string, '(A32)') processingrange%starttime
        case ( p_endtime )
          call get_string ( sub_rosa(subtree(2,son)), name_string )
          read(name_string, '(A32)') processingrange%endtime
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
          call forwardModelGlobalSetup ( son )
        case ( s_forwardModel )
          call decorate (son, AddForwardModelConfigToDatabase ( &
            & forwardModelConfigDatabase, ConstructForwardModelConfig ( son, vGrids ) ) )
        case ( s_vgrid )
          call decorate ( son, AddVGridToDatabase ( vGrids, &
            & CreateVGridFromMLSCFInfo ( name, son, l2gpDatabase ) ) )
        case ( s_l1brad )
          call l1bradSetup ( son, l1bInfo, F_FILE, &
          & MAXNUML1BRADIDS, ILLEGALL1BRADID )
        case ( s_l1boa )
          call l1boaSetup ( son, l1bInfo, F_FILE )

        case ( s_time )
          if ( timing ) then
            call sayTime
          else
            call cpu_time ( t1 )
            timing = .true.
          end if
        end select
      end if
    end do

   if( ECHO_GLOBAL_STNGS ) &
   & call dump_global_settings( l2pcf, processingRange, l1bInfo )

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

    ! Begin
    version = 1

    call output ( 'L1B database:', advance='yes' )
  
   if(associated(l1bInfo%L1BRADIDs)) then
    if ( num_l1b_files > 0 ) then
      do i = 1, num_l1b_files
!        returnStatus = Pgs_pc_getReference(l1bInfo%L1BRADIDs(i), version, &
!        & physicalFilename)
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

    call output ( 'L1OA file:', advance='yes' )
  
!    returnStatus = Pgs_pc_getReference(l1bInfo%L1BOAID, version, &
!      & physicalFilename)
!    if ( returnStatus == PGS_S_SUCCESS ) then
      call output ( 'fileid:   ' )
      call output ( l1bInfo%L1BOAID, advance='yes' )
      call output ( 'name:   ' )
      call output ( TRIM(l1bInfo%L1BOAFileName), advance='yes' )
!    else
!      call output ( '(file unknown)', advance='yes' )
!    end if

    call output ( 'Start Time:   ' )
    call output ( l2pcf%startutc, advance='yes' )

    call output ( 'End Time:   ' )
    call output ( l2pcf%endutc, advance='yes' )

    call output ( 'Start Time (tai):   ' )
    call output ( processingrange%starttime, advance='yes' )

    call output ( 'End Time:   ' )
    call output ( processingrange%endtime, advance='yes' )

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

  end subroutine dump_global_settings

end module GLOBAL_SETTINGS

! $Log$
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
