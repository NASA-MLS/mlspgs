module GLOBAL_SETTINGS

  use ForwardModelConfig, only: AddForwardModelConfigToDatabase, &
    & ForwardModelConfig_T
  use ForwardModelInterface, only: ConstructForwardModelConfig, &
    & ForwardModelGlobalSetup
  use HGrid, only: AddHGridToDatabase, CreateHGridFromMLSCFInfo, &
    & DestroyHGridDatabase, HGrid_T
  use INIT_TABLES_MODULE, only: L_TRUE, P_ALLOW_CLIMATOLOGY_OVERLOADS, &
    & P_INPUT_VERSION_STRING, P_OUTPUT_VERSION_STRING, P_VERSION_COMMENT, &
    & S_FORWARDMODEL, S_ForwardModelGlobal, S_TIME, S_VGRID
  use L2GPData, only: L2GPDATA_T
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Allocate
  use MoreTree, only: GET_FIELD_ID, GET_SPEC_ID
  use Output_m, only: Output
  use String_Table, only: Get_String
  use TOGGLES, only: GEN, LEVELS, SWITCHES, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATE, DECORATION, NODE_ID, NSONS, SUB_ROSA, SUBTREE
  use VGrid, only: CreateVGridFromMLSCFInfo, Dump
  use VGridsDatabase, only: AddVGridToDatabase, VGrid_T
  use TREE_TYPES, only: N_EQUAL, N_NAMED

  implicit NONE

  private

  public :: SET_GLOBAL_SETTINGS

  logical, public :: ALLOW_CLIMATOLOGY_OVERLOADS = .false.
  integer, public :: INPUT_VERSION_STRING = 0     ! Sub_rosa index
  integer, public :: OUTPUT_VERSION_STRING = 0    ! Sub_rosa index
  integer, public :: VERSION_COMMENT = 0          ! Sub_rosa index

!---------------------------- RCS Ident Info ---------------------------
  character (len=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!-----------------------------------------------------------------------

contains

  subroutine SET_GLOBAL_SETTINGS ( ROOT, ForwardModelConfigDatabase, &
    & VGrids, l2gpDatabase )

    integer, intent(in) :: ROOT    ! Index of N_CF node in abstract syntax tree
    type(ForwardModelConfig_T), dimension(:), pointer :: &
      & ForwardModelConfigDatabase
    type ( vGrid_T ), pointer, dimension(:) :: VGrids
    type ( l2gpData_T), dimension(:), pointer :: L2GPDATABASE
    
    integer :: I         ! Index of son of root
    integer :: NAME      ! Sub-rosa index of name of vGrid or hGrid
    integer :: SON       ! Son of root
    logical :: TIMING    ! For S_Time
    real :: T1, T2       ! For S_Time

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
      call output ( "Timing for MLSSignals = " )
      call output ( dble(t2 - t1), advance = 'yes' )
      timing = .false.
    end subroutine SayTime

  end subroutine SET_GLOBAL_SETTINGS

end module GLOBAL_SETTINGS

! $Log$
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
